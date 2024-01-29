#!/bin/bash 

# Test pbmm2 with PacBio long reads
#
# valeria
# srun -p medium -c 8 --time=7-00:00 -J safoPUVx_001-21_test_pbmm2 --mem=100G -o log/safoPUVx_001-21_test_pbmm2_%j.log /bin/sh 01_scripts/test_pbmm2.sh "safoPUVx_001-21" &

# VARIABLES
SAMPLE=$1
GENOME="03_genome/genome.fasta"
RAW_DATA_DIR="04_raw_data"
FILT_DIR="05_filtered"
CCS_BAM="$RAW_DATA_DIR/"$SAMPLE".ccs.bam" 
SUBREADS_BAM="$RAW_DATA_DIR/"$SAMPLE".subreads.bam" 

ALIGNED_DIR="06_aligned"

CPU=8

# 0. Create output dir
if [[ ! -d "$ALIGNED_DIR"/pbmm2 ]]
then
  mkdir $ALIGNED_DIR/pbmm2
fi

# 0. Index reference
#pbmm2 index $GENOME 03_genome/genome.fasta.mmi

# 1. Run pbmm2, sort and indem bam output
pbmm2 align $GENOME $CCS_BAM $ALIGNED_DIR/pbmm2/"$SAMPLE".ccs.bam --sort --sample $SAMPLE --bam-index BAI --preset CCS -j $(($CPU / 2)) -J $(($CPU / 2)) --log-level WARN  --rg "@RG\tID:$SAMPLE\tSM:$SAMPLE"



# try with subreads
#pbmm2 align $GENOME $SUBREADS_BAM $ALIGNED_DIR/pbmm2/"$SAMPLE".subreads.bam --sort --sample $SAMPLE --bam-index BAI --preset SUBREAD -j $(($CPU / 2)) -J $(($CPU / 2)) --log-level TRACE --rg "@RG\tID:$SAMPLE\tSM:$SAMPLE"

# try with trimmed subreads 
#pbmm2 align $GENOME QC/longQC/$SAMPLE/"$SAMPLE".subreads.trimmed.fastq $ALIGNED_DIR/pbmm2/"$SAMPLE".subreads.trimmed.bam --sort --sample $SAMPLE --bam-index BAI --preset SUBREAD -j $(($CPU / 2)) -J $(($CPU / 2)) --log-level TRACE --rg "@RG\tID:$SAMPLE\tSM:$SAMPLE"

