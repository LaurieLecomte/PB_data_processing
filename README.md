# ONT_data_processing

Pipeline for filtering and mapping Oxford Nanopore (ONT) reads

## Pipeline Overview

1. Concatenate raw fastq files for each sample : `00_cat_fastq_files.sh`
2. Plot summary statistics for each sample (optional) : `00_NanoPlot.sh` 
3. Filter reads according to lenght and minimum quality : `01_NanoFilt.sh` 
4. Map reads to the reference : `02_winnowmap.sh` 

## Prerequisites

### Files

* A reference genome named `genome.fasta` and its index (.fai) in `03_genome`
* Raw ccs bam files (or fastq) in `04_raw_data`, one file per sample.
* A sample IDs list (`02_infos/ind_PB.txt`), one ID per line

### Software

#### Required tools
* [SMRTLink tools](https://www.pacb.com/support/software-downloads/)
* [Winnowmap 2.03+](https://github.com/marbl/Winnowmap/releases/tag/v2.03)
* samtools 1.15
* GNU parallel

#### conda

1. Create env `PacBio_preprocess` from file : `conda create --name PacBio_preprocess --file PacBio_preprocess_env.txt`
2. Activate env : `conda activate PacBio_preprocess`

