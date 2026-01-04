import os
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from jenkspy import JenksNaturalBreaks
import pysam
import traceback
from tqdm import tqdm
import numpy as np
pysam_bam = None
import pandas as pd
import json
import argparse
from collections import Counter
import math
from math import comb
import io
import mysql.connector
from mysql.connector import errorcode


def load_csv_to_df(csv_tables_path, results_tables_list):
      print(f"Loading CSV tables from {csv_tables_path}...")
      csv_files = {f: os.path.join(csv_tables_path, f) for f in results_tables_list}
      
      missing_files = [f for f in csv_files if not os.path.isfile(f)]
      if missing_files:
          raise FileNotFoundError(f" CSV files ({str(missing_files)}) not found in directory {csv_tables_path}")
      else:
          print("All CSV files exist.")
    
      results_df_dict = dict()
      for name, path in csv_files.items():
          df = pd.read_csv(path)
          results_df_dict[name] = df
      
      return results_df_dict


def load_mysql_to_df(sql_db, results_tables_list):
      conn = sql_db.get("connection")
      if conn is None:
        raise ValueError("SQL database connection is not established.")
    
      results_df_dict = dict()
      for table in results_tables_list:
        try:
            query = f"SELECT * FROM `{table}`;"
            df = pd.read_sql(query, conn)
            results_df_dict[table] = df
        except Exception as e:
            raise RuntimeError(f"Failed to load table '{table}': {e}")

      return results_df_dict


def save_aggr_to_db(results_name, aggregated_df, sql_db):
    """
    Save aggregated DataFrame to MySQL table using mysql.connector.
    - Raises exception if table already exists
    - Automatically serializes dict/list columns to JSON
    """

    conn = sql_db["connection"]
    cursor = sql_db["cursor"]
    database = sql_db["database"]

    # Check if table exists
    cursor.execute(f"DROP TABLE IF EXISTS `{results_name}`")

    df = aggregated_df.replace({np.nan: None})

    # Detect JSON-like columns and serialize
    json_columns = []
    for col in df.columns:
        if df[col].apply(lambda x: isinstance(x, (dict, list))).any():
            df[col] = df[col].apply(
                lambda x: json.dumps(x) if isinstance(x, (dict, list)) else x
            )
            json_columns.append(col)

    # Generate CREATE TABLE statement
    columns_sql = []
    for col_name, dtype in zip(df.columns, df.dtypes):
        if col_name in json_columns:
            sql_type = "JSON"
        elif pd.api.types.is_integer_dtype(dtype):
            sql_type = "INT"
        elif pd.api.types.is_float_dtype(dtype):
            sql_type = "FLOAT"
        elif pd.api.types.is_bool_dtype(dtype):
            sql_type = "BOOLEAN"
        else:
            sql_type = "VARCHAR(255)"
        columns_sql.append(f"`{col_name}` {sql_type}")

    create_table_sql = f"""
        CREATE TABLE `{results_name}` (
            {', '.join(columns_sql)}
        );
    """

    cursor.execute(create_table_sql)
    conn.commit()

    # Prepare INSERT
    placeholders = ", ".join(["%s"] * len(df.columns))
    insert_sql = f"""
        INSERT INTO `{results_name}`
        ({', '.join(f'`{c}`' for c in df.columns)})
        VALUES ({placeholders})
    """

    data_tuples = [tuple(row) for row in df.itertuples(index=False, name=None)]

    # Define a chunk size
    chunk_size = 1000
    print(f"Inserting data in chunks of {chunk_size}...")
    
    for i in range(0, len(data_tuples), chunk_size):
        chunk = data_tuples[i : i + chunk_size]
        try:
            cursor.executemany(insert_sql, chunk)
            conn.commit() # Commit after each chunk to keep memory usage low
        except mysql.connector.Error as err:
            print(f"Error inserting chunk starting at row {i}: {err}")
            conn.rollback()
            raise

    print(
        f"Data saved to MySQL table '{results_name}' "
        f"({len(data_tuples)} rows total)."
    )


def build_haplotypes_dict(lengths, threshold):
    """
    Clusters haplotype lengths by relative difference threshold,
    handling zero-length values without division errors.
    Returns dict mapping cluster-mean -> {'c': count, 'p': frequency}.
    """
    def cluster_floats(values, threshold):
        if not values:
            return []
    
        floats = sorted(float(x) for x in values)
        clusters = []
        current_cluster = [floats[0]]
        current_min = floats[0]
        
        for candidate in floats[1:]:
            # Special handling if current_min is zero
            if current_min == 0:
                if candidate == 0:
                    current_cluster.append(candidate)
                    continue
                else:
                    clusters.append(current_cluster)
                    current_cluster = [candidate]
                    current_min = candidate
                    continue
            # Normal relative threshold check
            if (candidate - current_min) / current_min < threshold:
                current_cluster.append(candidate)
            else:
                clusters.append(current_cluster)
                current_cluster = [candidate]
                current_min = candidate
        clusters.append(current_cluster)
        return clusters
    
    # Flatten all haplotype lengths    
    clusters = cluster_floats(lengths, threshold)
    result = {}
    for cluster in clusters:
        count = len(cluster)
        mean_len = sum(cluster) / count
        result[mean_len] = {'c': count, 'p': count / len(lengths)}
    return result
         

def rarefy_K_from_dict(hap_dict, rarefied_n):
    """
    Rarefies number of distinct haplotypes to sample size rarefied_n
    from a dict of haplotype counts.
    """
    hap_counts = [v['c'] for v in hap_dict.values()]
    total_n = sum(hap_counts)
    if total_n < rarefied_n:
        return None  # Not enough data to rarefy

    rarefied_K = 0
    for count in hap_counts:
        # Probability that haplotype is NOT picked at all
        p_not_drawn = comb(total_n - count, rarefied_n) / comb(total_n, rarefied_n)
        rarefied_K += (1 - p_not_drawn)
    return rarefied_K
        

def calculate_intron_variability(args):
    """
    intron_id, gene_id, table_for_stats, threshold_for_len_difference, rarefied_n
    Aggregate intron length and read counts from all samples for a given intron_id.
    """
    (intron_id, gene_id, table_for_stats, rarefied_n, threshold_for_len_difference) = args
    
    try:
        samples = list()
        length_haplotypes_list = list()
        sample_reads = dict()
        sample_haplotypes = dict()
        het = list()
        x_haplotypes = None
        x_length_difference = None
        x_fc = None
        
        for sample_name, sample_df in shared_df_dict.items():
            if intron_id not in sample_df.index:
                continue
                
            row = sample_df.loc[intron_id]
            
            if row['trusted'] != 1:
                continue
            
            if table_for_stats == sample_name:
                x_haplotypes = [row['length_mean1'], row['length_mean2']]
                x_length_difference = max(x_haplotypes)-min(x_haplotypes)
                x_fc = max(x_haplotypes)/min(x_haplotypes)
                            
            if row['length_mean1'] != row['length_mean2']:
                het.append(sample_name)
            samples.append(sample_name)
            length_haplotypes_list.append(row['length_mean1'])
            length_haplotypes_list.append(row['length_mean2'])
            sample_reads[sample_name] = row['intron_lengths']
            sample_haplotypes[sample_name] = [row['length_mean1'], row['length_mean2']]
            
                       
        if len(samples) == 0:
            return None
            
        n = len(samples)*2
        clustered_haplotypes = build_haplotypes_dict(length_haplotypes_list, threshold_for_len_difference)
        clustered_haplotypes_list = [key for key in clustered_haplotypes.keys()]
        k = len(clustered_haplotypes_list)
        krare = rarefy_K_from_dict(clustered_haplotypes, rarefied_n)
        ho = len(het)/len(samples)

        data = {
            "id_intron": intron_id,
            "Analyzed_samples": samples,
            "Intron_lengths_for_individual_reads": sample_reads,
            "Length_haplotypes_for_individual_samples": sample_haplotypes,
            "Length_haplotypes_list": length_haplotypes_list,
            "Clustered_length_haplotypes": clustered_haplotypes,
            "Clustered_length_haplotypes_list": clustered_haplotypes_list,
            "Number_of_haplotypes_n": n ,
            "Haplotype_richness_K": k,
            "Rarefied_n6_haplotype_richness_K": krare,
            "Observed_heterozygosity_Ho": ho,
            "SRX14125033_haplotypes": x_haplotypes,
            f"SRX14125033_haplotype_length_difference": x_length_difference,
            f"SRX14125033_haplotype_length_fold_change": x_fc
        }
        
        return data
        
    except Exception as e:
        raise RuntimeError(f"Error processing intron {intron_id} (gene {gene_id})") from e


def aggregate_intron_results_table(results_df_dict, results_name, genes_dict, table_for_stats, rarefied_n, threshold_for_len_difference=0.10, proc=8):
    """
    Aggregate intron length results across multiple samples using genes_dict for annotation.
    """
    print(f"### Aggregating intron results for {results_name} ###")

    # Index each results table by IntronID for fast lookup
    indexed_results = dict()
    for t_name, df in results_df_dict.items():
        if isinstance(df, pd.DataFrame):
            indexed_results[t_name] = df.set_index('id_intron').sort_index()
        else:
            raise ValueError("Each value in results_df_dict must be a pandas DataFrame")
            
    # Aggregate over all introns in genes_dict keys
    tasks = [
        (intron_id, gene_id, table_for_stats, rarefied_n, threshold_for_len_difference) 
        for gene_id in genes_dict 
        for intron_id in genes_dict[gene_id]["introns"]
    ]
    chunksize = max(1, 10 * proc)
    aggregated_results = []
    print("Aggregating intron data in parallel...")
    
    with ProcessPoolExecutor(max_workers=proc, initializer=init_worker2, initargs=(genes_dict, indexed_results),) as executor:
        results = executor.map(calculate_intron_variability, tasks, chunksize=chunksize,)    

        for res in tqdm(results, total=len(tasks)):
            if res is not None:
                aggregated_results.append(res)

    print(f"Aggregation completed for {len(aggregated_results)} introns.")
    return aggregated_results
    
    
def init_worker(bam_path, genes_dict, static_params):
    global pysam_bam, shared_genes, w_params
    pysam_bam = pysam.AlignmentFile(bam_path, "rb")
    shared_genes = genes_dict
    w_params = static_params
    
def init_worker2(gene_dict, df_dict):
    global shared_gene_dict
    global shared_df_dict
    shared_gene_dict = gene_dict 
    shared_df_dict = df_dict 

    
def extract_read_ids_from_bam(chrom, start, end):
    # global pysam_bam
    read_ids = set()
    for read in pysam_bam.fetch(chrom, start, end):
        read_id = read.query_name
        if ";" not in read_id:
            read_ids.add(read.query_name)
    return read_ids
    

def blast_exons_combined(exon1_seq, exon1_id, exon2_seq, exon2_id, subject_ids, long_read_blastdb, threads):
    """
    Runs a single BLAST process for both exons by streaming to stdin.
    """
    # Prepare FASTA format string for both exons
    fasta_input = f">{exon1_id}\n{exon1_seq}\n>{exon2_id}\n{exon2_seq}\n"
    
    # BLAST+ requires -seqidlist to be a physical file
    with tempfile.NamedTemporaryFile(mode="w+", delete=True) as id_file:
        id_file.write("\n".join(subject_ids) + "\n")
        id_file.flush()

        cmd = [
            "blastn", "-task", "blastn", 
            "-db", long_read_blastdb, 
            "-seqidlist", id_file.name,
            "-outfmt", "6 qseqid qlen qstart qend sseqid slen sstart send pident length bitscore qseq sseq",
            "-soft_masking", "false", 
            "-max_target_seqs", "10000", 
            "-num_threads", str(threads), 
            "-dust", "no"
        ]

        try:
            process = subprocess.Popen(
                cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE, text=True
            )
            stdout, stderr = process.communicate(input=fasta_input)
            
            if process.returncode != 0:
                raise Exception(f"BLAST error: {stderr}")
            
            return stdout
        except Exception as e:
            raise Exception(f"Error in blast_exons_combined: {e}")
            
    

def get_jenks_means(lst, diff_thr):
    if len(set(lst)) > 1:
        jnb = JenksNaturalBreaks(2)
        jnb.fit(lst)
        mean1 = np.average(jnb.groups_[0])
        mean2 = np.average(jnb.groups_[1])

        if abs(mean2 - mean1) < (diff_thr * mean2):
            num = len(lst)
            mean = sum(lst) / len(lst)
            sd = (sum([((x - mean) ** 2) for x in lst]) / num) ** 0.5
            final = [mean, mean, num, num, sd, sd]
        else:
            num1 = len(jnb.groups_[0])
            num2 = len(jnb.groups_[1])
            sd1 = (sum([((x - mean1) ** 2) for x in jnb.groups_[0]]) / num1) ** 0.5
            sd2 = (sum([((x - mean2) ** 2) for x in jnb.groups_[1]]) / num2) ** 0.5
            final = [mean1, mean2, num1, num2, sd1, sd2]
    else:
        num = len(lst)
        mean = sum(lst) / len(lst)
        sd = (sum([((x - mean) ** 2) for x in lst]) / num) ** 0.5
        final = [mean, mean, num, num, sd, sd]

    return final
    
    
def intron_len_calculation(exon1_blastOut, exon2_blastOut, bit_score, diff_thr, margin):
    length_computed = False
    
    if (exon1_blastOut == "") or (exon2_blastOut == ""):
        selected_sseqids = {}
    else:
        # blast output to dictionary
        column_names = ["qseqid", "qlen", "qstart", "qend", "sseqid", "slen", "sstart", "send", "pident", "length", "bitscore", "qseq", "sseq"]
        columns_to_extract = ["qseqid", "qlen", "qstart", "qend", "slen", "sstart", "send", "pident", "length", "bitscore", "qseq", "sseq"]
        column_indices = {col: i for i, col in enumerate(column_names)}
        sseqid_index = column_indices["sseqid"]
        exon1_blast_dic = {}
        
        rows = exon1_blastOut.strip().split('\n')
        for row in rows:
            if not row:
                continue
            values = row.split('\t')
            if len(values) < len(column_names):
                continue

            sseqid = values[sseqid_index]
            value_dict = {col: values[column_indices[col]] for col in columns_to_extract}
            exon1_blast_dic[sseqid] = value_dict
        exon2_blast_dic = {}
        for row in exon2_blastOut.strip().split('\n'):
            values = row.split('\t')
            sseqid = values[sseqid_index]
            value_dict = {col: values[column_indices[col]] for col in columns_to_extract}
            exon2_blast_dic[sseqid] = value_dict
    
        # extracting reads with both exons blast hits
        selected_sseqids = dict()
        for sseqid in exon1_blast_dic:
            if sseqid in exon2_blast_dic:
                record1 = exon1_blast_dic[sseqid]
                record2 = exon2_blast_dic[sseqid]
                
                # Check bitscore > bit_score and length within qlen +/- margin for both records
                if (float(record1['bitscore']) > bit_score and float(record2['bitscore']) > bit_score and
                    (int(record1['qlen']) - margin) < int(record1['length']) < (int(record1['qlen']) + margin) and
                    (int(record2['qlen']) - margin) < int(record2['length']) < (int(record2['qlen']) + margin)): 
                    if (int(record1['sstart']) < int(record1['send']) and 
                        int(record2['sstart']) < int(record2['send']) and 
                        int(record1['send']) < int(record2['sstart'])
                       ):
                        i_len = int(record2['sstart']) - int(record1['send']) - 1
                        i_e1 = record1['sstart'] + "-" + record1['send']
                        i_e2 = record2['sstart'] + "-" + record2['send']
                    elif (int(record1['sstart']) > int(record1['send']) and 
                          int(record2['sstart']) > int(record2['send']) and 
                          int(record2['sstart']) < int(record1['send'])
                         ):
                        i_len = int(record1['send']) - int(record2['sstart']) - 1
                        i_e1 = record1['sstart'] + "-" + record1['send'] + "r"
                        i_e2 = record2['sstart'] + "-" + record2['send'] + "r"
                    else:
                        continue
                    selected_sseqids[sseqid] = [int(i_len), sseqid + "#e1:" + i_e1 + "#e2:" + i_e2]     

    # store results
    if len(selected_sseqids) > 0:
        length_computed = True
        read_IDs = [value[1] for value in selected_sseqids.values()]
        lengths = [value[0] for value in selected_sseqids.values()]
        jenks_means = get_jenks_means(lengths, diff_thr)
        num_reads = len(read_IDs)
        length_mean1 = jenks_means[0]
        length_mean2 = jenks_means[1]
        num_reads1 = jenks_means[2]
        num_reads2 = jenks_means[3]
        length_sd1 = jenks_means[4]
        length_sd2 = jenks_means[5]
    else:
        read_IDs = None
        lengths = None
        num_reads = 0
        length_mean1 = None
        length_mean2 = None
        num_reads1 = 0
        num_reads2 = 0
        length_sd1 = None
        length_sd2 = None
    return (length_computed,
            read_IDs, lengths, num_reads, length_mean1, length_mean2, num_reads1, num_reads2, length_sd1, length_sd2)
            

def is_not_nonzero_number(x):
    return not (isinstance(x, (int, float)) and not math.isnan(x) and x != 0)


def check_single_triple(read_lens_list, len1_val, len2_val, diff, num_reads1_val, num_reads2_val, intron_len_thr, read_num_thr):
    """
    True if all read lengths are within at least one of the intervals:
       [len1_val*(1 - diff), len1_val*(1 + diff)] or
       [len2_val*(1 - diff), len2_val*(1 + diff)].
       AND difference between read numbers is les than 10x
       AND max length > intron_len_thr
       AND number of reads per intron > read_num_thr
    If any read length is outside *both* intervals, return False (untrust).
    """

    # Check if len1_val and len2_val are nonzero numbers
    if is_not_nonzero_number(len1_val) or is_not_nonzero_number(len2_val):
        return False

    # check if all reads are close to allele lengths
    len1_low  = len1_val * (1 - diff)
    len1_high = len1_val * (1 + diff)
    len2_low  = len2_val * (1 - diff)
    len2_high = len2_val * (1 + diff)

    for rl in read_lens_list:
        outside_len1 = (rl < len1_low) or (rl > len1_high)
        outside_len2 = (rl < len2_low) or (rl > len2_high)

        # If it's outside BOTH intervals => untrust
        if outside_len1 and outside_len2:
            return False  # untrust
    
    unbalanced = min(num_reads1_val, num_reads2_val)*10 < max(num_reads1_val, num_reads2_val)
    if unbalanced:
        return False  # untrust
    
    # check if sufficient number of reads
    if len(read_lens_list) < read_num_thr:
        return False  # untrust

    # check if sufficient intron length 
    if max(read_lens_list) < intron_len_thr:
        return False  # untrust        
    
    return True  # otherwise trust
        
            
def get_intron_length_statistics(gene_id):
    """
    Processes a gene's introns by blasting both exons against long reads in one go.
    """
    gene = shared_genes[gene_id]
    w = w_params  # Shorthand for config
    
    if 'introns' not in gene or not gene['introns']:
        return False
        
    out_rows_list = []
    
    for intron_name, intron in gene['introns'].items():
        read_ids = extract_read_ids_from_bam(intron['chrom'], intron['start'], intron['end'])
        
        if len(read_ids) <= 1:
            exon1_blastOut = exon2_blastOut = ""
        else:
            # 1. Perform a single BLAST for both exons
            raw_blast_out = blast_exons_combined(
                intron['exon1_seq'], intron['exon1_id'],
                intron['exon2_seq'], intron['exon2_id'],
                read_ids, w["nanopore_db_path"], 1
            )
            
            # 2. Parse results to keep ONLY the best hit per Query-Subject pair
            # This replaces the 'sort -u' shell logic safely for multiple queries
            best_hits = {}
            for line in raw_blast_out.strip().split('\n'):
                if not line: continue
                cols = line.split('\t')
                q_id, s_id, length, bitscore = cols[0], cols[4], int(cols[9]), float(cols[10])
                
                key = (q_id, s_id)
                if key not in best_hits:
                    best_hits[key] = line
                else:
                    # Logic: Higher bitscore wins; tie-break with alignment length
                    cur_bit = float(best_hits[key].split('\t')[10])
                    cur_len = int(best_hits[key].split('\t')[9])
                    if bitscore > cur_bit or (bitscore == cur_bit and length > cur_len):
                        best_hits[key] = line

            # 3. Separate the hits by exon ID
            e1_lines = [line for (q, s), line in best_hits.items() if q == intron['exon1_id']]
            e2_lines = [line for (q, s), line in best_hits.items() if q == intron['exon2_id']]
            
            exon1_blastOut = "\n".join(e1_lines)
            exon2_blastOut = "\n".join(e2_lines)

        # 4. Length calculation and filtering
        (length_computed, read_IDs, lengths, num_reads, l_mean1, l_mean2, 
         n_reads1, n_reads2, l_sd1, l_sd2) = intron_len_calculation(
            exon1_blastOut, exon2_blastOut, w["bit_score"], w["diff_thr"], w["margin"]
        )

        trusted = 1 if check_single_triple(
            lengths, l_mean1, l_mean2, w["diff_thr"], n_reads1, n_reads2, 
            w["min_intron_len"], w["min_read_num"]
        ) else 0

        out_rows_list.append({
            "IntronID": intron_name, "Position": intron_name, "Gene_name": gene["gene_name"],
            "Species": w["species_name"], "Gene_ID": gene_id, "exon1_name": intron['exon1_id'],
            "exon2_name": intron['exon2_id'], "intron_lengths": lengths, "read_ids": read_IDs,
            "num_reads": num_reads, "length_mean1": l_mean1, "length_mean2": l_mean2,
            "num_reads1": n_reads1, "num_reads2": n_reads2, "length_sd1": l_sd1,
            "length_sd2": l_sd2, "trusted": trusted
        })
                          
    return out_rows_list


def filter_by_blacklist(gene_res_list, intron_blacklist):
    """
    Filter gene_res_list by removing entries whose intron ID
    is present in intron_blacklist.

    Parameters
    ----------
    gene_res_list : list of dict
        Output list from IntronLength.work()
    intron_blacklist : list or set
        List of intron IDs to exclude

    Returns
    -------
    list of dict
        Filtered gene_res_list
    """
    blacklist_set = set(intron_blacklist)

    filtered = [
        res for res in gene_res_list
        if res["IntronID"] not in blacklist_set
    ]
    
    print(f"{len(filtered)} remained from {len(gene_res_list)} - {(len(filtered)/len(gene_res_list))*100}%")
    return filtered


def work(species_name: str, sample_name: str, nanopore_fasta: str, genome_fasta_file: str, genes_dict: dict, bit_score: int = 50, margin: int = 20, diff_thr: float = 0.1, gen_lim = None, gene_list = None, proc = 8, intron_blacklist = None, max_read_num_diff = 10, min_intron_len = 100, min_read_num = 4):
    """
    """
    print(f"### calculating varible introns for species: {species_name}, sample: {sample_name} ###")
    nanopore_fasta_no_suffix = os.path.splitext(nanopore_fasta)[0]
    minimap_out_filename = nanopore_fasta_no_suffix + "_" + species_name + ".minimap.bam"
    nanopore_db_path = nanopore_fasta_no_suffix
    
    # creating blast database
    blastdb_extensions = [".ndb"]
    if all(os.path.exists(nanopore_db_path + ext) for ext in blastdb_extensions):
        print(f"BLAST database already exists: {nanopore_db_path}")
    else:
        cmd = ["makeblastdb", "-in", nanopore_fasta, "-dbtype", "nucl", "-parse_seqids", "-out", nanopore_db_path]
        print(f"running command: {str(cmd)}")
        subprocess.run(cmd, check=True)
        print("BLAST database created:", nanopore_db_path)
        
    # mapping reads to reference
    if os.path.exists(minimap_out_filename):
        print(f"minimap results already present in: {minimap_out_filename}")
    else:
        print(f"starting minimap for {nanopore_fasta}")
        cmd = (
            f"minimap2 --secondary=no -a -t 18 -x map-ont {genome_fasta_file} {nanopore_fasta} "
            f"| samtools view -Sb - "
            f"| samtools sort -o {minimap_out_filename} - "
            f"&& samtools index {minimap_out_filename}")
        print(f"running command: {str(cmd)}")
        result = subprocess.run(cmd, shell=True, text=True, check=True)
        print("minimap bam created:", minimap_out_filename)
    
    # Process pool for intron length calculation per genes
    print("intron length calculation (gene looping)")
    gene_res_list = []
    excs = []

    # Prepare the list of genes to process
    if gene_list:
        genes_to_process = {gid: gene for gid, gene in genes_dict.items() if gid in gene_list}
    else:
        genes_to_process = genes_dict

    # Apply gen_lim if specified
    if gen_lim is not None:
        genes_to_process = dict(list(genes_to_process.items())[:gen_lim])

    print(f"Submitting {len(genes_to_process)} genes for intron length calculation...")

    # Submit jobs
    static_params = dict(
        species_name=species_name,
        minimap_out_filename=minimap_out_filename,
        nanopore_db_path=nanopore_db_path,
        bit_score=bit_score,
        diff_thr=diff_thr,
        margin=margin,
        max_read_num_diff=max_read_num_diff,
        min_intron_len=min_intron_len,
        min_read_num=min_read_num,
    )
    with ProcessPoolExecutor(max_workers=proc, initializer=init_worker, initargs=(minimap_out_filename, genes_dict, static_params)) as pool:
        futures = {pool.submit(get_intron_length_statistics, gene_id): gene_id for gene_id in genes_to_process}

        # Collect results as they finish
        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing genes"):
            gid = futures[future]
            try:
                res = future.result()
                if res:
                    gene_res_list.extend(res)
            except Exception as e:
                excs.append(traceback.format_exception(None, e, e.__traceback__))
                tqdm.write(f"Exception in gene {gid}: {e}")

    print(f"{len(excs)} exceptions and {len(gene_res_list)} results from {len(genes_to_process)} genes")
    print("all results retrived")            
    
    if intron_blacklist:
        gene_res_list = filter_by_blacklist(gene_res_list, intron_blacklist)
    
    return excs, gene_res_list
            
