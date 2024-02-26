# PB_data_processing

Pipeline for mapping PacBio CCS reads to a reference genome.

## Pipeline Overview

1. Pre-compute kmers: `01_compute_kmers.sh`
2. Convert ccs bam files to fasta, then map reads to the reference for each sample: `02_winnowmap.sh` 


### Additional scripts

Other scripts targeting a specific step or operation conducted in one of the main scripts or allowing additional analyses are provided in the `01_scripts/utils` subdirectory.

* `01_scripts/utils/generate_ccs.sh` : Generate a .ccs file from .subreads file for a given sample.
* `01_scripts/utils/longQC_ccs.sh` : Get quality metrics from a .ccs file for a given sample.

Older scripts used for development or debugging purposes are stored in the `01_scripts/archive` folder for future reference if needed. These are not meant to be used in their current state and may be obsolete.



## Prerequisites

### Files

* A reference genome named `genome.fasta` and its index (.fai) in `03_genome`
* Raw ccs bam files (or fastq) in `04_raw_data`, one file per sample.
* A sample IDs list (`02_infos/ind_PB.txt`), one ID per line

### Software

#### Required tools
* [SMRTLink tools](https://www.pacb.com/support/software-downloads/)
* [Winnowmap 2.03+](https://github.com/marbl/Winnowmap/releases/tag/v2.03)
* GNU parallel


#### For Manitou users
Custom conda environments are required for running `smrt-tools` and `winnowmap`, as these programs are not available on Manitou; See the [Conda environment preparation](#conda-environment-preparation) section below. 

#### For users working with other computing clusters and servers
The program versions specified in this pipeline refer to the versions available on IBIS' bioinformatics servers when this pipeline was built in 2021-2022, and are likely not available on all other servers. 
Please add a '#' at the beginning of each line in the `#LOAD REQUIRED MODULES` section in each script (or remove these lines), and follow the [Conda environment preparation](#conda-environment-preparation) to create custom conda environments with correct program versions and dependencies.
A R installation is also required.


## Detailed Walkthrough

For running each script, copy the `srun` command from the script's header to the terminal and adjust parameters (memory, partition, time limit) if necessary.  
The header also features a brief description of the script's contents. 


### Conda environment preparation

#### SV calling environments (`winnowmap` + `smrt-tools`)
From the main directory, run `conda create --name PacBio_preprocess --file PacBio_preprocess_env.txt`

### Main pipeline

#### 1. Pre-compute kmers (`01_compute_kmers.sh`)

#### 2. Map CCS reads to the reference genome (`02_winnowmap.sh`)


