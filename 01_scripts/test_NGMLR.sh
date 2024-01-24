#!/bin/bash 

# Test NGMLR with PacBio long reads
#
# valeria
# srun -p medium -c 8 --time=7-00:00 -J safoPUVx_001-21_test_NGMLR --mem=100G -o log/safoPUVx_001-21_test_NGMLR_%j.log /bin/sh 01_scripts/test_NGMLR.sh "safoPUVx_001-21" &

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
if [[ ! -d "$ALIGNED_DIR"/NGMLR ]]
then
  mkdir $ALIGNED_DIR/NGMLR
fi

# 0. Index reference
#pbmm2 index $GENOME 03_genome/genome.fasta.mmi

# CCS
ngmlr -t $CPU -r $GENOME -q $RAW_DATA_DIR/"$SAMPLE".ccs.fasta.gz -o $ALIGNED_DIR/NGMLR/"$SAMPLE".ccs.bam -x pacbio --rg-id "$SAMPLE" --rg-sm "$SAMPLE"



