# Intron Variability

This repository contains a Python-based framework for investigating genomic structural variation, specifically focusing on **Intron Length Polymorphism** using long-read sequencing data (Oxford Nanopore Technologies). 

The code serves as a methodological component for the scientific publication in *Genome Biology and Evolution*.

## 1. General Description
The pipeline automates the extraction of genomic features from standard annotation formats, simulates long-read data for quality control, and performs high-throughput BLAST-based screening to measure intron lengths across cohorts.

## 2. Reference
If you use this code in your research, please cite:
> *[Full Paper Citation to be inserted: Author(s). (Year). Title. Genome Biology and Evolution, Volume(Issue), pages.]*

---

## 3. Prerequisites & Installation

### System Dependencies
The following tools must be installed and accessible in your system `PATH`:
* `minimap2`
* `blastn` (BLAST+)
* `samtools`
* `MySQL Server`

### Environment Setup
A Conda `environment.yml` file is provided in the root directory to manage Python dependencies. To create and activate the environment, run:

```bash
conda env create -f environment.yml
conda activate genome_pipeline
```

The environment includes essential libraries such as pysam, pybedtools, gffutils, mysql-connector-python, pandas, bioframe, and jenkspy

## 4. Pipeline Methods

The workflow is managed via `execution.py` and configured through a YAML file. Key modules include:

### **Genome Initialization & Pre-processing**
* **Annotation Parsing**: Converts GTF/GFF3 files into a standardized Python dictionary. It handles gene merging for overlapping features and filters out micro-exons (default <10bp).
* **Coordinate Validation**: Cross-references the annotation against the FASTA index to ensure all coordinates are within valid contig bounds.

### **Intron Length Calculation**
* **Read Mapping**: Uses `minimap2` to align Nanopore reads to the reference genome.
* **Fragment Analysis**: Identifies reads containing hits for both flanking exons. The distance between the end of Exon 1 and the start of Exon 2 on the long read determines the intron length.
* **Allele Calling**: Employs **Jenks Natural Breaks** optimization to cluster read lengths into distinct alleles, accounting for sequencing error profiles.

### **Aggregation & Statistics**
* **Population Metrics**: Aggregates results across a cohort to calculate:
    * **Haplotype Richness ($K$)**: Total distinct clustered alleles.
    * **Rarefied $K$**: Richness normalized to a specific sample size ($n=6$) to allow comparison between samples.
    * **Observed Heterozygosity ($H_o$)**: Proportion of analyzed samples showing multi-allelic status.

---

## 5. Output Files and Data Structure

### **Flat Files (CSV & BED)**

| File Suffix | Description |
| :--- | :--- |
| `_exons.bed` | 6-column BED file of raw exons extracted from the initial annotation. |
| `_merged_and_filtered_exons.bed` | Final BED file after merging overlapping genes and filtering short exons. |
| `_intron_lengths.csv` | Per-sample results containing measured introns, read IDs, and calculated means. |
| `_intron_lengths_aggr.csv` | Cohort-wide summary with calculated $K$, $H_o$, and rarefaction data. |

### **MySQL Database Tables**
The pipeline organizes data into a relational schema for efficient querying:

* **`{species}_chromosomes`**: Metadata for contigs, including length and `tr_coverage` (repeat fraction).
* **`{species}_genes`**: Table of merged gene boundaries, orientations, and names.
* **`{species}_exons`**: Coordinates of individual exons linked back to parent genes.
* **`{species}_introns`**: Target intron coordinates, flanking exons, and `tr_coverage`.
* **`{species}_{sample_name}`**: Sample-specific results including JSON-encoded read lengths and trusted status.
* **`{results_name}`**: The final aggregated table containing population statistics across the cohort.

---

## 6. Usage

1. **Configure**: Edit the `example.yaml` file to define your file paths, SQL credentials, and analysis thresholds.
2. **Run**:

```bash
python execution.py your_config.yaml
```
