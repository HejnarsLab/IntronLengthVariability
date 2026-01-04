import sys
import yaml
from src import genome_class as gc


def load_config(config_path):
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


if __name__ == "__main__":

    # ---------------------
    # CLI ARGUMENT
    # ---------------------
    if len(sys.argv) != 2:
        print("Usage: python execution.py <config.yaml>")
        sys.exit(1)

    cfg = load_config(sys.argv[1])
    steps = cfg["steps"]

    # ---------------------
    # GENERAL
    # ---------------------
    proc = cfg["proc"]
    gen_lim = cfg["gen_lim"]
    gene_list = cfg["gene_list"]

    # ---------------------
    # GENOME INIT (always required)
    # ---------------------
    mySpecies = gc.GenomeSample(
        species_name=cfg["species_name"],
        outfile_prefix=cfg["outfile_prefix"],
        fasta_path=cfg["fasta_path"],
        gtf_path=cfg["gtf_path"],
        gtf_db_path=cfg["gtf_db_path"],
        bed_file_path=cfg["bed_file_path"],
        annot_format=cfg["annot_format"],
        min_exon_len=cfg["min_exon_len"]
    )

    mySpecies.get_summary()

    # ---------------------
    # ADD SQL DATABASE
    # ---------------------
    if steps.get("add_sql_database", False):
        sql = cfg["sql_database"]
        mySpecies.add_sql_database(
            host=sql["host"],
            user=sql["user"],
            password=sql["password"],
            database=sql["database"],
            port=sql["port"],
        )

    # ---------------------
    # CREATE SQL DATABASE
    # ---------------------
    if steps.get("create_sql_database", False):
        mySpecies.create_sql_database()
    
            
    # ---------------------
    # CREATE INTRON BLACKLIST
    # ---------------------
    if steps.get("create_intron_blacklist", False):
        blacklist = cfg["create_intron_blacklist"]
        calc = cfg["calculate_intron_lengths"]

        mySpecies.create_intron_blacklist(
            nanopore_fasta=blacklist["sim_nanopore_fasta"],
            coverage=blacklist["coverage"],
            min_len=blacklist["min_len"],
            max_len=blacklist["max_len"],
            seed=blacklist["seed"],
            max_read_num_diff=calc["max_read_num_diff"],
            min_intron_len=calc["min_intron_len"],
            min_read_num=calc["min_read_num"],
            bit_score=calc["bit_score"],
            margin=calc["margin"],
            diff_thr=calc["diff_thr"],
            proc=proc,
            gen_lim=gen_lim,
            gene_list=gene_list,
            results_table_path=blacklist["results_table_path"]
        )


    # ---------------------
    # CALCULATE INTRON LENGTHS
    # ---------------------
    if steps.get("calculate_intron_lengths", False):
        calc = cfg["calculate_intron_lengths"]

        for sample_name, nanopore_fasta in zip(
            calc["sample_name_list"],
            calc["nanopore_fasta_list"],
        ):
            mySpecies.calculate_intron_lengths(
                sample_name,
                nanopore_fasta,
                store_to_db=calc["store_to_db"],
                max_read_num_diff=calc["max_read_num_diff"],
                min_intron_len=calc["min_intron_len"],
                min_read_num=calc["min_read_num"],
                bit_score=calc["bit_score"],
                margin=calc["margin"],
                diff_thr=calc["diff_thr"],
                proc=proc,
                gen_lim=gen_lim,
                gene_list=gene_list,
            )

    # ---------------------
    # AGGREGATE RESULTS
    # ---------------------
    if steps.get("aggregate_intron_results", False):
        agg = cfg["aggregate_intron_results"]
        calc = cfg["calculate_intron_lengths"]

        mySpecies.aggregate_intron_results(
            results_tables_list=agg["results_tables_list"],
            results_name=agg["results_name"],
            table_for_stats=agg["table_for_stats"],
            threshold_for_len_difference=calc["diff_thr"],
            proc=proc,
            store_to_db=agg["store_to_db"],
            csv_tables_path=agg["csv_tables_path"],
            rarefied_n=agg["rarefied_n"],
        )

    mySpecies.get_summary()
