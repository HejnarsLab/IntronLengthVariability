####
# testovani zatim u tvorby blacklistu
# funkce generuje same nulove hodnoty: excs, gene_res_list = IntronLength.work(self.species_name, "simulated", nanopore_fasta, self.fasta_path, self.genes, bit_score = bit_score, margin = margin, diff_thr = diff_thr, gen_lim = gen_lim, gene_list = gene_list, proc = proc, intron_blacklist = None, max_read_num_diff = max_read_num_diff, min_intron_len = min_intron_len, min_read_num = min_read_num)
####

# libraries
import os
import re
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import pysam
import gffutils
import pybedtools
import mysql.connector
from mysql.connector import errorcode
import json
import subprocess
import csv
import pandas as pd
import random

from . import intron_lengths as IntronLength
from . import exon_stuttering as ExonStuttering

# functions
def validate_paths(file_path):
    """Checks if the provided file path is valid"""
    if not os.path.exists(file_path):
        raise Exception(f"file not found at: {file_path}")


def gff_db_to_dict(gff_db, annot_format):
    """
    Convert gffutils object to gene dictionary. gtf and gff modes enabled.
    Note: GTF mode uses coordinate containment (limit) to find exons for transcripts.
    """
    gene_dict = dict()
    annot_format = annot_format.lower() # Ensure format is lower case
    
    # Determine features to iterate over based on format
    if annot_format == "gtf":
        gff_genes = gff_db.all_features(featuretype="transcript")
        # gff_gene_count = gff_db.count_features_of_type(featuretype="transcript")
    elif annot_format == "gff":
        gff_genes = gff_db.all_features(featuretype="gene")
        # gff_gene_count = gff_db.count_features_of_type(featuretype="gene")
    else:
        raise ValueError(f"annot_format parameter must be 'gtf' or 'gff', not '{annot_format}'")

    # Iterate and extract exon data
    for gene in gff_genes:
        try:
            # Handle attribute access differences:
            if annot_format == "gtf":
                gene_id = gene.attributes.get("gene_id", ["UNKNOWN_ID"])[0].replace('|', '_').replace(',', '_')
                gene_name = gene.attributes.get("gene_name", [gene_id])[0].replace('|', '_').replace(',', '_')
            else: # gff
                gene_id = gene.attributes.get("ID", ["UNKNOWN_ID"])[0].replace('|', '_').replace(',', '_')
                gene_name = gene.attributes.get("Name", [gene_id])[0].replace('|', '_').replace(',', '_')

            
            # GTF/GFF Exon Retrieval logic
            if annot_format == "gtf":  # Using coordinate containment fot GTF (not idela)
                childrens = gff_db.all_features(
                    featuretype="CDS", 
                    limit=(gene.chrom, gene.start, gene.stop), 
                    strand=gene.strand, 
                    order_by='start', 
                    completely_within=True
                )
            else: # Use standard, reliable children() method for GFF
#                childrens = list()
#                for mrna in gff_db.children(gene, featuretype="mRNA", order_by="start"):
#                    cds_list = list(gff_db.children(mrna, featuretype="CDS", order_by="start"))
#                    childrens.extend(cds_list)
                childrens = gff_db.children(
                    gene, 
                    featuretype='CDS', 
                    order_by='start')
            
            # Retriving exon coordinates
            exon_starts = list()
            exon_ends = list()
            exon_chroms = list()
            
            for i in childrens:
                exon_starts.append(int(i.start)-1)
                exon_ends.append(int(i.stop))
                exon_chroms.append(i.chrom)
                
                # test if start is not higher than stop
                if int(i.start)>int(i.stop):
                    raise Exception(f"{i.start} not >= {i.stop} in {gene_id} exon")
            
            # Check if exon coordinates are present
            if not (len(exon_starts) > 0 and len(exon_starts) == len(exon_ends) == len(exon_chroms)):
                continue

            # Check if all exons are on a single contig
            if len(set(exon_chroms)) != 1: 
                raise Exception(f"Gene/Transcript is located on more than 1 contig: {set(exon_chroms)}. ID: {gene_id}")
            
            # gene coords
            gene_pos = [exon_chroms[0], min(exon_starts), max(exon_ends)]
            # gene_pos = [gene.chrom, min([int(gene.start), int(gene.stop)])-1, max([int(gene.start), int(gene.stop)])]
            #if gene_pos[0] != exon_chroms[0] or gene_pos[1] != min(exon_starts) or gene_pos[2] != max(exon_ends):
            #    raise Exception(f"Gene/Transcript coordinates do not correspond to exon coords. Gene: {gene_id}, Pos: {gene_pos}, chr:{exon_chroms}, starts: {exon_starts}, ends: {exon_ends}")
            
            gene_dict[gene_id] = {
               "gene_name": gene_name,
               "gene_position": gene_pos,
               "exon_chroms": exon_chroms, 
               "exon_starts": exon_starts,
               "exon_ends": exon_ends,
               "strand": str(gene.strand)
            }
            
        except Exception as e:
            raise Exception(f"Something wrong during gtf/gff parsing in {gene.id if hasattr(gene, 'id') else 'N/A'}). error: {e}")

    # Final check
    if len(gene_dict) == 0:
        raise Exception(f"Empty result of gff_db_to_dict function. No valid features found in '{annot_format}' format.")
    
    return gene_dict
            

def validate_annotations_against_fasta(genes_dict, fasta_contigs):
    """
    Checks if all exon coordinates in the genes dictionary are valid against the FASTA index.
    """
    print("\n--- Starting Genome Coordinate Validation ---")
    for gene_id, data in genes_dict.items():
        gene_name = data["gene_name"]
        
        # Check 1: Gene/Exon chromosome must exist in FASTA index
        # Since all exons should be on the same chrom (checked in gff_db_to_dict), we check the first one
        exon_chrom = data["exon_chroms"][0] 
        if exon_chrom not in fasta_contigs:
            raise Exception(
                f"Validation Error: Chromosome '{exon_chrom}' for gene {gene_id} (Name: {gene_name}) "
                f"is not found in the FASTA index.")
        contig_length = fasta_contigs[exon_chrom]
        
        # Check 2: All exon coordinates must be within the contig length
        for i, (start, end) in enumerate(zip(data["exon_starts"], data["exon_ends"])):
            # BED coordinates are 0-based, start must be >= 0. End must be <= contig_length.
            
            if start < 0:
                 raise Exception(
                    f"Validation Error: Exon {i+1} of gene {gene_id} has invalid start coordinate ({start}). "
                    f"Coordinate must be >= 0."
                )
            
            if end > contig_length:
                raise Exception(
                    f"Validation Error: Exon {i} of gene {gene_id} has invalid end coordinate ({end}). "
                    f"Contig '{exon_chrom}' length is {contig_length}, but exon end is {end}."
                )
    
    print("Genome Coordinate Validation successful: All annotations are within FASTA contig bounds.")
    
    
def save_exons_to_bed(genes_dict, output_filename):
    """
    Writes exon coordinates from the genes dictionary to a 6-col BED file.
    Format: chrom\tstart\tend\tgene_name|gene_id\tscore\tstrand
    Note: We use 0-based start for true BED compatibility.
    """
    with open(output_filename, 'w') as outfile:
        # Iterate over all gene IDs and their data
        for gene_id, data in genes_dict.items():
            gene_name = data["gene_name"]
            gene_strand = data["strand"]
            
            # Iterate through the exons
            for chrom, start, end in zip(data["exon_chroms"], data["exon_starts"], data["exon_ends"]):
                custom_name = f"{gene_name}|{gene_id}"
                score = '0'
                line = f"{chrom}\t{start}\t{end}\t{custom_name}\t{score}\t{gene_strand}\n"
                outfile.write(line)
                

def filter_exons_by_size(gene_dict, genomeIndex, min_exon_len = 10):
    """Must be at least 10 bp long (end - start)."""
    print("filtering exons <10bp and adding introns dictionary to gene dict.")
    gene_filtered_dict = dict()
    for gene_id in tqdm(gene_dict.keys(), total = len(gene_dict)):
        exon_chroms = list()
        exon_starts = list()
        exon_ends = list()
        exon_ids = list()
        for chrom, start, end in zip(gene_dict[gene_id]["exon_chroms"], gene_dict[gene_id]["exon_starts"], gene_dict[gene_id]["exon_ends"]):
            if (int(end) - int(start)) >= min_exon_len:
                exon_chroms.append(chrom)
                exon_starts.append(int(start))
                exon_ends.append(int(end))
                exon_ids.append(f"{str(chrom)}:{str(start)}-{str(end)}")
        
        # test if exons are present
        if len(exon_chroms) == 0:
            continue
        
        # reset gene coordinate
        gene_coord = [exon_chroms[0], min(exon_starts), max(exon_ends)]

        # create intron dict
        introns = dict()
        for i, (chrom, start, end, exid) in enumerate(zip(exon_chroms, exon_starts, exon_ends, exon_ids)):
            # create exon seq
            seq = str(genomeIndex.fetch(chrom, start, end))      
            
            # first exon
            if i == 0:
                start_1 = start
                end_1 = end
                seq_1 = seq
                id_1 = exid 
                continue
            
            # defining coordinates
            int_start = end_1
            int_end = start
            ex1_start = start_1
            ex1_end = end_1
            ex2_start = start
            ex2_end = end
            intron_id = f"{str(chrom)}:{str(int_start)}-{str(int_end)}"
            ex1_id = id_1
            ex2_id = exid
            
            # exon seqs
            exon_seq = seq_1
            exon_seq2 = seq
            
            # set for other itteration
            start_1 = start
            end_1 = end
            seq_1 = seq
            id_1 = exid 

            introns[intron_id] = {
                "chrom": chrom,
                "start": int_start,
                "end": int_end,
                "exon1_start": ex1_start,
                "exon1_end": ex1_end,
                "exon2_start": ex2_start,
                "exon2_end": ex2_end,
                "exon1_id": ex1_id,
                "exon2_id": ex2_id,
                "exon1_seq": exon_seq,
                "exon2_seq": exon_seq2}

        # insert to gene_filtered_dict
        gene_filtered_dict[gene_id] = {
            "gene_name": gene_dict[gene_id]["gene_name"],
            "gene_position": gene_coord,
            "exon_chroms": exon_chroms, 
            "exon_starts": exon_starts,
            "exon_ends": exon_ends,
            "exon_ids": exon_ids,
            "strand": gene_dict[gene_id]["strand"], 
            "introns": introns}        
    
    print("exons filtered")    
    return gene_filtered_dict
    
    
def merge_genes(gene_dict, exon_bed_path):
    """merge overlapping genes to one gene entry"""
    print(f"merging overlapping genes - {len(gene_dict)} genes")
    gene_merged_dict = dict()
    
    # gene coord bed object creation and merging
    gene_coord_generator = (pybedtools.Interval(val['gene_position'][0], val['gene_position'][1], val['gene_position'][2], name=key, score="0", strand=val['strand'], otherfields=[val['gene_name']]) for key, val in gene_dict.items())
    gene_bed = pybedtools.BedTool(gene_coord_generator).sort().saveas()
    merged_gene_bed = gene_bed.merge(s=False, c='4,5,6,7', o='distinct,sum,distinct,distinct', d=-1).saveas()

    # exon coord bed object creation and merging
    exons_bed = pybedtools.BedTool(exon_bed_path).sort().saveas()
    merged_exons_bed = exons_bed.merge(s=False, d=-1).saveas()
    
    # recreating gene_dict
    exons_genes_intersect_bed = merged_exons_bed.intersect(merged_gene_bed, wa=True, wb=True).sort().saveas()
    for feature in tqdm(exons_genes_intersect_bed):
        chrom, start, end, chrom_g, start_g, end_g, gene_id_com, score, strand_com, gene_name_com  = feature.fields
        strands = strand_com.split(',')
        if len(set(strands)) > 1 and len(set(strands)) != 0:
            strand = '.'
        else:
            strand = str(strands[0])
        try:
            gene_merged_dict[gene_id_com]["exon_chroms"].append(chrom)
            gene_merged_dict[gene_id_com]["exon_starts"].append(int(start))
            gene_merged_dict[gene_id_com]["exon_ends"].append(int(end))
        except KeyError:
            gene_merged_dict[gene_id_com] = {
               "gene_name": gene_name_com,
               "gene_position": [chrom_g, int(start_g), int(end_g)],
               "exon_chroms": [chrom], 
               "exon_starts": [int(start)],
               "exon_ends": [int(end)],
               "strand": strand}                

    print(f"merging overlapping genes finished - {len(gene_dict)-len(gene_merged_dict)} genes merged")
    return gene_merged_dict
      

def bed_to_gene_dict(bed_path):
    """Converts a bed file into the gene dictionary format"""
    genes_dict = {}
    gene_bed = pybedtools.BedTool(bed_path).sort().saveas()
    
    for feature in tqdm(gene_bed):
        chrom, start, end, name, score, strand = feature.fields

        try:
            gene_id, gene_name = name.split('|')
        except ValueError:
            gene_id, gene_name = [name, name]

        try:
            genes_dict[gene_id]["exon_chroms"].append(chrom)
            genes_dict[gene_id]["exon_starts"].append(int(start))
            genes_dict[gene_id]["exon_ends"].append(int(end))
            genes_dict[gene_id]["gene_position"] = [chrom, min(genes_dict[gene_id]["exon_starts"]), max(genes_dict[gene_id]["exon_ends"])]
        except KeyError:
            genes_dict[gene_id] = {
               "gene_name": gene_name,
               "gene_position": [chrom, int(start), int(end)],
               "exon_chroms": [chrom], 
               "exon_starts": [int(start)],
               "exon_ends": [int(end)],
               "strand": str(strand)}   

    return genes_dict


def create_species_sql_schema_from_connection(species_name: str, sequence_dict: dict, genes: dict, sql_db: dict, force: bool = False):
    """
    Create and populate species-specific tables using an existing SQL connection.
    """

    T_CHR = f"{species_name}_chromosomes"
    T_GEN = f"{species_name}_genes"
    T_EX  = f"{species_name}_exons"
    T_INT = f"{species_name}_introns"

    conn = sql_db["connection"]
    cursor = sql_db["cursor"]
    database = sql_db["database"]

    def table_exists(name):
        cursor.execute("""
            SELECT COUNT(*) FROM information_schema.tables
            WHERE table_schema=%s AND table_name=%s
        """, (database, name))
        return cursor.fetchone()[0] > 0

    tables = [T_INT, T_EX, T_GEN, T_CHR]

    if not force:
        for t in tables:
            if table_exists(t):
                raise Exception(
                    f"Table '{t}' already exists. Use force=True."
                )
    else:
        for t in tables:
            cursor.execute(f"DROP TABLE IF EXISTS {t}")
            
    print("creating tables")
    # Create tables
    cursor.execute(f"""
        CREATE TABLE {T_CHR} (
            id_chr VARCHAR(255) PRIMARY KEY,
            name VARCHAR(255),
            length BIGINT,
            tr_coverage FLOAT,
            properties_dict JSON,
            filter VARCHAR(255),
            class VARCHAR(255),
            status VARCHAR(255),
            note TEXT
        )
    """)

    cursor.execute(f"""
        CREATE TABLE {T_GEN} (
            id_gene VARCHAR(500) PRIMARY KEY,
            chr VARCHAR(255),
            start BIGINT,
            end BIGINT,
            orientation VARCHAR(255),
            name VARCHAR(500),
            filter VARCHAR(255),
            class VARCHAR(255),
            status VARCHAR(255),
            note TEXT,
            FOREIGN KEY (chr) REFERENCES {T_CHR}(id_chr)
        )
    """)

    cursor.execute(f"""
        CREATE TABLE {T_EX} (
            id_exon VARCHAR(255) PRIMARY KEY,
            chr VARCHAR(255),
            start BIGINT,
            end BIGINT,
            id_gene VARCHAR(500),
            filter VARCHAR(255),
            class VARCHAR(255),
            status VARCHAR(255),
            note TEXT,
            FOREIGN KEY (chr) REFERENCES {T_CHR}(id_chr),
            FOREIGN KEY (id_gene) REFERENCES {T_GEN}(id_gene)
        )
    """)

    cursor.execute(f"""
        CREATE TABLE {T_INT} (
            id_intron VARCHAR(255) PRIMARY KEY,
            chr VARCHAR(255),
            start BIGINT,
            end BIGINT,
            upstream_exon VARCHAR(255),
            downstream_exon VARCHAR(255),
            id_gene VARCHAR(500),
            length INT,
            tr_coverage FLOAT,
            filter VARCHAR(255),
            class VARCHAR(255),
            status VARCHAR(255),
            note TEXT,
            FOREIGN KEY (chr) REFERENCES {T_CHR}(id_chr),
            FOREIGN KEY (upstream_exon) REFERENCES {T_EX}(id_exon),
            FOREIGN KEY (downstream_exon) REFERENCES {T_EX}(id_exon),
            FOREIGN KEY (id_gene) REFERENCES {T_GEN}(id_gene)
        )
    """)
    print("tables created")

    # 1. Prepare data lists for batch insertion
    chr_data = [(chrom, chrom, length) for chrom, length in sequence_dict.items()]
    
    gene_data = []
    exon_data = []
    intron_data = []

    print("Preparing data for batch insertion...")
    for gene_id, gene in tqdm(genes.items()):
        chrom, start, end = gene["gene_position"]
        gene_data.append((gene_id, chrom, start, end, gene["strand"], gene["gene_name"]))

        for ex_s, ex_e, ex_id in zip(gene["exon_starts"], gene["exon_ends"], gene["exon_ids"]):
            exon_data.append((ex_id, chrom, ex_s, ex_e, gene_id))

        for intron_id, intron in gene["introns"].items():
            length = intron["end"] - intron["start"] + 1
            intron_data.append((intron_id, intron["chrom"], intron["start"], intron["end"], 
                               intron["exon1_id"], intron["exon2_id"], gene_id, length))

    # 2. Execute Batch Inserts
    try:
        # --- DIAGNOSTIC CHECK ---
        # Check if any gene_id in your list is too long before sending to MySQL
        oversized_genes = [(item[0], len(item[0])) for item in gene_data if len(str(item[0])) > 500]
        
        if oversized_genes:
            print(f"\n[!] Found {len(oversized_genes)} gene IDs that are too long (> 500 chars):")
            for long_id, length in oversized_genes[:5]: # Print first 5
                print(f"Length: {length} | ID: {long_id}")
            print("... (Stopping execution to prevent SQL error)")
            return # Exit the function early
        # ------------------------
        print(f"Inserting {len(chr_data)} chromosomes...")
        cursor.executemany(f"INSERT INTO {T_CHR} (id_chr, name, length) VALUES (%s,%s,%s)", chr_data)
        
        print(f"Inserting {len(gene_data)} genes...")
        cursor.executemany(f"INSERT INTO {T_GEN} (id_gene, chr, start, end, orientation, name) VALUES (%s,%s,%s,%s,%s,%s)", gene_data)
        
        print(f"Inserting {len(exon_data)} exons...")
        cursor.executemany(f"INSERT INTO {T_EX} (id_exon, chr, start, end, id_gene) VALUES (%s,%s,%s,%s,%s)", exon_data)
        
        print(f"Inserting {len(intron_data)} introns...")
        cursor.executemany(f"INSERT INTO {T_INT} (id_intron, chr, start, end, upstream_exon, downstream_exon, id_gene, length) VALUES (%s,%s,%s,%s,%s,%s,%s,%s)", intron_data)
        
        conn.commit()
        print("Done! All data committed.")
    except mysql.connector.Error as err:
        print(f"Error: {err}")
        conn.rollback() # Undo changes if one fails



def simulate_nanopore_reads(pysam_fasta_index, fasta_dict, output_fasta: str, coverage: float, min_len: int, max_len: int, seed: int = 42):
    """
    Simulate Nanopore-like reads by sampling genome subsequences using seqkit.
    """
    random.seed(seed)
    contigs = list(fasta_dict.keys())
    lengths = list(fasta_dict.values())
    total_len = sum(lengths)
    # Length-weighted contig sampling
    probs = [l / total_len for l in lengths]
    
    if min_len >= max_len:
        raise ValueError("min_len must be smaller than max_len")

    avg_read_len = (min_len + max_len) / 2
    number_reads = int((coverage * total_len) / avg_read_len)

    if number_reads < 1:
        raise ValueError("Calculated number of reads < 1; adjust coverage or read lengths")
    
    BASES = "ACGT"
    
    def reverse_complement(seq: str) -> str:
        return seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]
    
    def substitute_N(seq: str) -> str:
        return "".join(random.choice(BASES) if b in "Nn" else b for b in seq)

    with open(output_fasta, "w") as out:
        read_id = 1
        for read_id in tqdm(range(1, number_reads + 1), desc="Simulating reads"):
            while True:
                contig = random.choices(contigs, weights=probs, k=1)[0]
                contig_len = fasta_dict[contig]

                read_len = random.randint(min_len, max_len)
                if read_len > contig_len:
                    continue

                start = random.randint(0, contig_len - read_len)

                seq = pysam_fasta_index.fetch(contig, start, start + read_len)

                # Replace Ns
                seq = substitute_N(seq)

                # Random strand
                if random.random() < 0.5:
                    seq = reverse_complement(seq)

                out.write(
                    f">read_{read_id}_{contig}:{start}-{start + read_len}\n"
                )
                out.write(seq + "\n")
                break


def load_blacklist_from_csv(csv_path):
    """
    Loads the simulated results CSV and reconstructs the intron_blacklist.
    """
    print(f"Loading existing results from {csv_path}...")
    
    # Load the CSV into a DataFrame
    df = pd.read_csv(csv_path)

    # check variable
    mask = (df['length_mean1'] != df['length_mean2']) & (df['length_mean1'].notnull()) & (df['length_mean1']>0)
    invalid_count = df['length_mean1'].isnull().sum() + (df['length_mean1'] == 0).sum()
    blacklist = df.loc[mask, 'id_intron'].tolist()
    
    print(f"Blacklist reconstructed: {len(blacklist)} introns identified as variable, {invalid_count} null/zero skipped.")
    return blacklist
          


# classes
class GenomeSample:
    """
    Represents a single genome assembly with its sequence, annotation,
    and optional analysis results like repeat content etc.
    
    Params: outfile_prefix, fasta_path [gtf_path, annot_format | gtf_db_path, annot_format | bed_file_path]
        species_name - name of species
        outfile_prefix - path with file prefix for output files
        fasta_path - Path to the FASTA file containing the genomic sequence.
        gtf_path - Path to the annotation file (GTF or GFF).
        gtf_db_path - Path to stored gtf database file
        bed_file_path - Path to bed file with exon coordinates and 4th column with gene names in format: "gene_name|gene_id"
        annot_format - Specifies the format of the annotation file ('gtf' or 'gff').
    
    Attributes:
        self.species_name 
        self.fasta_path
        self.gtf_path 
        self.gtf_db_path
        self.annot_format 
        self.input_bed_path 
        self.bed_path - path to unmerged exon bed file
        self.filtered_bed_path - path to merged and filtered exon bed file 
        self.outfile_prefix 
        self.sql_db - must be created by method 
        self.simulated_nanopore_fasta_path - must be created by method 
        self.intron_blacklist - must be created by method 
        self.fasta_index - pysam.FastaFile object
        self.sequence_dict
        self.genes - dictionary of gene annotations where overlapped genes were merged and <10b exons were excluded
        
    
    Methods:
        add_sql_database(host: str, user: str, password: str, database: str, port: int = 3306)
        create_sql_database(force: bool = False)
        simple_repeat_coverage(bed_file_path)
        exon_stuttering_screen(blacklist_bed_file: str | None, hit_length_threshold: int, min_intron_len: int, store_to_db, bitscore_thr: float, min_copy_num: int, proc=8, store_shuffled=False)
        create_intron_blacklist(nanopore_fasta: str | None = None, coverage: float = 100.0, min_len: int = 5000, max_len: int = 50000, seed: int = 42)
        calculate_intron_lengths(sample_name, nanopore_fasta, store_to_db = True, max_read_num_diff = 10, min_intron_len = 100, min_read_num = 4, bit_score = 50, margin = 20, diff_thr = 0.1, gen_lim = None, gene_list = None, proc = 8)
        aggregate_intron_results(results_tables_list, results_name, threshold_for_len_difference=0.10, proc=8, store_to_db=True, csv_tables_path=False, table_for_stats, rarefied_n = 6)
        calculate_exon_stuttering()
        calculate_TR_coverage()
        get_summary()
    
    Out files:
        gff database file (if bed file not provided)
        gene bed file
    """

    def __init__(self, species_name: str, outfile_prefix: str, fasta_path: str, gtf_path: str | None = None, gtf_db_path: str | None = None, annot_format: str | None = None, bed_file_path: str | None = None, min_exon_len: int = 10):
        """init code"""
        print("class building")
        # Store metadata
        self.species_name = species_name.replace(" ", "_")
        self.fasta_path = fasta_path
        self.gtf_path = gtf_path
        self.gtf_db_path = gtf_db_path
        self.annot_format = annot_format.lower()
        self.input_bed_path = bed_file_path
        self.bed_path = outfile_prefix + "_" + species_name + "-exons.bed"
        self.filtered_bed_path = outfile_prefix + "_" + species_name + "-merged_and_filtered_exons.bed"
        self.outfile_prefix = outfile_prefix
        self.sql_db = None
        self.simulated_nanopore_fasta_path = None
        self.intron_blacklist = None
        self.fasta_index = None
        self.sequence_dict = None
        self.genes = None
        
        
        # Input Path Validation
        validate_paths(fasta_path)
        validate_paths(gtf_path)

        # 1. Load fasta index to attribute
        validate_paths(fasta_path)
        try:
            # fasta index object
            print(f"creating fasta index from {fasta_path}")
            self.fasta_index = pysam.FastaFile(fasta_path)
            
            # Pre-fetch chromosome lengths
            self.sequence_dict = dict(zip(self.fasta_index.references, self.fasta_index.lengths))
            
            # print results
            conts50 = list()
            for i, name in enumerate(self.sequence_dict.keys()):
                if i < 50:  
                    conts50.append(name)
                else:
                    break
            print(f"\nFirst 50 contig names: {conts50}\n\nFASTA file successfully indexed. Total records: {len(self.fasta_index.references)}.")
            
        except Exception as e:
            raise Exception(f"Failed to create fasta index for '{fasta_path}': {e}")

        # 2. Creat gene annotations
        # 2.1 create gene dict
        if bed_file_path is not None:
            genes_raw = bed_to_gene_dict(bed_file_path)
            print(f"gene dictionary created. {len(genes_raw)} genes \n")
        
        elif gtf_db_path is not None:
            validate_paths(gtf_db_path)
            gtf_db = gffutils.interface.FeatureDB(dbfn=gtf_db_path)
            print(f"gff database loaded\n")
            genes_raw = gff_db_to_dict(gtf_db, self.annot_format)
            print(f"gene dictionary created. {len(genes_raw)} genes \n")
            
        elif gtf_path is not None:
            validate_paths(gtf_path)
            print(f"gff database creating from {gtf_path}")
            gtf_db_path = gtf_path + ".db"
            if not os.path.exists(gtf_db_path):
                gffutils.create_db(gtf_path, dbfn=gtf_db_path, force_gff=True, verbose=False, merge_strategy="create_unique")
                print(f"gff database created: {gtf_db_path}\n")
            else:
                print(f"gff database already present in: {gtf_db_path}\n")
            gtf_db = gffutils.interface.FeatureDB(dbfn=gtf_db_path)
            print(f"gff database loaded\n")                           
            genes_raw = gff_db_to_dict(gtf_db, self.annot_format)
            print(f"gene dictionary created. {len(genes_raw)} genes \n")
            
        else:
            raise Exception(f"no gtf_db_path, gtf_path, or bed_file_path provided")

        validate_annotations_against_fasta(genes_raw, self.sequence_dict) # validate gene dict
        save_exons_to_bed(genes_raw, self.bed_path) # safe to bed file
        print(f"gene dictionary saved to {self.bed_path}\n")

        # 2.2 create merged and filtered gene dict
        merged_genes = merge_genes(genes_raw, self.bed_path)
        self.genes = filter_exons_by_size(merged_genes, self.fasta_index, min_exon_len = min_exon_len)
        save_exons_to_bed(self.genes, self.filtered_bed_path) # safe to bed file
        print(f"gene dictionary saved to {self.filtered_bed_path}\n")


    def add_sql_database(self, host: str, user: str, password: str, database: str, port: int = 3306):
        """
        Attach to an existing MySQL database.
        """
        conn = mysql.connector.connect(host=host, user=user, password=password, database=database, port=port)
        cursor = conn.cursor()

        self.sql_db = {"connection": conn, "cursor": cursor, "database": database}

        print("SQL database connection established.")



    def create_sql_database(self, force: bool = False):
        """
        Create species-specific SQL schema using an existing DB connection.
        """
        if self.sql_db is None:
            raise Exception("SQL database not attached. Call add_sql_database() first.")

        create_species_sql_schema_from_connection(
            species_name=self.species_name,
            sequence_dict=self.sequence_dict,
            genes=self.genes,
            sql_db=self.sql_db,
            force=force
        )
    
    
    def create_intron_blacklist(self, nanopore_fasta: str | None = None, coverage: float = 100.0, min_len: int = 5000, max_len: int = 50000, seed: int = 42, max_read_num_diff = 10, min_intron_len = 100, min_read_num = 4, bit_score = 50, margin = 20, diff_thr = 0.1, gen_lim = None, gene_list = None, proc = 8, results_table_path: str | None = None):
        """
        Create blacklist of introns where length_mean1 != length_mean2.
        """
        if results_table_path:
            self.intron_blacklist = load_blacklist_from_csv(results_table_path)
        else:
            # 1. Simulate reads if FASTA not provided
            if nanopore_fasta is None:
                print("No nanopore FASTA provided, simulating reads")
                nanopore_fasta = (f"{self.outfile_prefix}_{self.species_name}-simulated.fa")

                simulate_nanopore_reads(self.fasta_index, self.sequence_dict, output_fasta=nanopore_fasta, coverage=coverage, min_len=min_len, max_len=max_len, seed=seed)
            else:
                print(f"Using provided nanopore FASTA: {nanopore_fasta}")

            self.simulated_nanopore_fasta_path = nanopore_fasta

            # 2. Run intron length calculation
            excs, gene_res_list = IntronLength.work(self.species_name, "simulated", nanopore_fasta, self.fasta_path, self.genes, bit_score = bit_score, margin = margin, diff_thr = diff_thr, gen_lim = gen_lim, gene_list = gene_list, proc = proc, intron_blacklist = None, max_read_num_diff = max_read_num_diff, min_intron_len = min_intron_len, min_read_num = min_read_num)

            # 3. Build blacklist
            blacklist = [res["IntronID"] for res in gene_res_list if (res["length_mean1"] != res["length_mean2"] and not IntronLength.is_not_nonzero_number(res["length_mean1"]))]
            self.intron_blacklist = blacklist
        
            # CSV writing
            outfile_path = self.outfile_prefix + "_" + self.species_name + "-simulated_intron_lengths.csv"
            
            with open(outfile_path, "w", newline="") as csvfile:
                fieldnames = ["id_intron", "species", "gene_name", "id_gene", "upstream_exon", "downstream_exon", "intron_lengths", "read_ids", "num_reads", "length_mean1", "length_mean2", "num_reads1", "num_reads2", "length_sd1", "length_sd2", "trusted"]
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()
                
                for res in gene_res_list:
                    writer.writerow({
                        "id_intron": res["IntronID"],
                        "species": res["Species"],
                        "gene_name": res["Gene_name"],
                        "id_gene": res["Gene_ID"],
                        "upstream_exon": res["exon1_name"],
                        "downstream_exon": res["exon2_name"],
                        "intron_lengths": json.dumps(res["intron_lengths"]),
                        "read_ids": json.dumps(res["read_ids"]),
                        "num_reads": res["num_reads"],
                        "length_mean1": res["length_mean1"],
                        "length_mean2": res["length_mean2"],
                        "num_reads1": res["num_reads1"],
                        "num_reads2": res["num_reads2"],
                        "length_sd1": res["length_sd1"],
                        "length_sd2": res["length_sd2"],
                        "trusted": res["trusted"]
                    })

            print(f"Intron length results stored in table {outfile_path}")
    
    
    def calculate_intron_lengths(self, sample_name, nanopore_fasta, store_to_db = True, max_read_num_diff = 10, min_intron_len = 100, min_read_num = 4, bit_score = 50, margin = 20, diff_thr = 0.1, gen_lim = None, gene_list = None, proc = 8):
        if store_to_db and (self.sql_db is None or "connection" not in self.sql_db):
            raise Exception("SQL database not initialized. Run create_sql_database().")
        if self.intron_blacklist is None:
            raise Exception("intron_blacklist not computed. Run create_intron_blacklist().")
        
        if store_to_db:    
            conn = self.sql_db["connection"]
            cursor = self.sql_db["cursor"]
            results_table = f"{self.species_name}_{sample_name}"
            
            cursor.execute(f"DROP TABLE IF EXISTS `{results_table}`")

            cursor.execute(f"""
                CREATE TABLE IF NOT EXISTS `{results_table}` (
                    id_intron VARCHAR(255) PRIMARY KEY,
                    intron_lengths JSON,
                    read_ids JSON,
                    num_reads INT,
                    length_mean1 FLOAT,
                    length_mean2 FLOAT,
                    num_reads1 INT,
                    num_reads2 INT,
                    length_sd1 FLOAT,
                    length_sd2 FLOAT,
                    trusted INT,
                    FOREIGN KEY (id_intron) REFERENCES {self.species_name}_introns(id_intron)
                )
            """)
        
        print(f"intron length calculation for {sample_name}")
        outfile_path = self.outfile_prefix + "_" + sample_name + "-intron_lengths.csv"
        excs, gene_res_list = IntronLength.work(self.species_name, sample_name, nanopore_fasta, self.fasta_path, self.genes, bit_score = bit_score, margin = margin, diff_thr = diff_thr, gen_lim = gen_lim, gene_list = gene_list, proc = proc, intron_blacklist = self.intron_blacklist, max_read_num_diff = max_read_num_diff, min_intron_len = min_intron_len, min_read_num = min_read_num)
        
        # Open CSV file for writing
        with open(outfile_path, "w", newline="") as csvfile:
            fieldnames = ["id_intron", "species", "gene_name", "id_gene", "upstream_exon", "downstream_exon", "intron_lengths", "read_ids", "num_reads", "length_mean1", "length_mean2", "num_reads1", "num_reads2", "length_sd1", "length_sd2", "trusted"]
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            
            for res in gene_res_list:
                writer.writerow({
                    "id_intron": res["IntronID"],
                    "species": res["Species"],
                    "gene_name": res["Gene_name"],
                    "id_gene": res["Gene_ID"],
                    "upstream_exon": res["exon1_name"],
                    "downstream_exon": res["exon2_name"],
                    "intron_lengths": json.dumps(res["intron_lengths"]),
                    "read_ids": json.dumps(res["read_ids"]),
                    "num_reads": res["num_reads"],
                    "length_mean1": res["length_mean1"],
                    "length_mean2": res["length_mean2"],
                    "num_reads1": res["num_reads1"],
                    "num_reads2": res["num_reads2"],
                    "length_sd1": res["length_sd1"],
                    "length_sd2": res["length_sd2"],
                    "trusted": res["trusted"]
                })
            
                if store_to_db:
                    cursor.execute(f"INSERT INTO `{results_table}` VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)",
                        (res["IntronID"],
                        json.dumps(res["intron_lengths"]),
                        json.dumps(res["read_ids"]),
                        res["num_reads"],
                        res["length_mean1"],
                        res["length_mean2"],
                        res["num_reads1"],
                        res["num_reads2"],
                        res["length_sd1"],
                        res["length_sd2"],
                        res["trusted"]))
            conn.commit()
        print(f"Intron length results stored in table `{results_table}`.")
        


    def aggregate_intron_results(self, results_tables_list, results_name, table_for_stats, threshold_for_len_difference=0.10, proc=8, store_to_db=True, csv_tables_path=False, rarefied_n = 6):
        """
        Wrapper function to aggregate intron results and store them as CSV or DB tables.
        """
        # Load tables from CSV if path is given
        if csv_tables_path:
            results_df_dict = IntronLength.load_csv_to_df(csv_tables_path, results_tables_list)
        else:
            if self.sql_db is None or "connection" not in self.sql_db:
                raise Exception("SQL database not initialized. Run create_sql_database().")            
            results_df_dict = IntronLength.load_mysql_to_df(self.sql_db, results_tables_list)

        # Run core aggregation function
        aggregated_results = IntronLength.aggregate_intron_results_table(results_df_dict, results_name, self.genes, table_for_stats, rarefied_n, threshold_for_len_difference=threshold_for_len_difference, proc=proc)

        # Convert aggregated results to DataFrame for saving
        aggregated_df = pd.DataFrame(aggregated_results)

        # Save CSV
        csv_output_path = self.outfile_prefix + "_" + results_name + "-intron_lengths_aggr.csv"
        aggregated_df.to_csv(csv_output_path, index=False)
        print(f"Aggregated results saved to {csv_output_path}")

        # Optionally store to DB
        if store_to_db:
            if self.sql_db is None or "connection" not in self.sql_db:
                raise Exception("SQL database not initialized. Run create_sql_database().")
            IntronLength.save_aggr_to_db(results_name, aggregated_df, self.sql_db)
            print(f"Aggregated results saved to DB table: {results_name}")
            
         
         
    def calculate_exon_stuttering(self):
        """ """
        print("need to be done")
        
    
    def calculate_TR_coverage(self):
        """ """
        print("need to be done")
 
             
    def get_summary(self) -> dict:
        """Returns a summary of the available data for the sample."""
        return {
            "species_name": self.species_name,
            "fasta_path": self.fasta_path,
            "gtf_path": self.gtf_path,
            "gtf_db_path": self.gtf_db_path,
            "annot_format": self.annot_format,
            "input_bed_path": self.input_bed_path,
            "bed_path": self.bed_path,
            "filtered_bed_path": self.filtered_bed_path,
            "outfile_prefix": self.outfile_prefix,
            "sql_db": self.sql_db,
            "simulated_nanopore_fasta_path": self.simulated_nanopore_fasta_path,
            "number of introns in blacklist": None if self.intron_blacklist is None else len(self.intron_blacklist),
            "number of contigs": len(self.sequence_dict),
            "number of genes": len(self.genes),
            "number of exons": len([intron_id for gene_id in self.genes for intron_id in self.genes[gene_id]["introns"]]),
            "example of one key from self.genes": next(iter(self.genes.keys())),
            "example of one item from self.genes": next(iter(self.genes.values()))}
