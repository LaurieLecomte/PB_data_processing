#!/bin/bash 

# Test pbmm2 with PacBio long reads
#
# valeria
# srun -p ibis_large -c 6 --time=21-00:00 -J test_pbmm2 --mem=70G -o log/test_pbmm2_%j.log /bin/sh 01_scripts/test_pbmm2.sh "safoCHIm_2359-21" &

# VARIABLES
SAMPLE=$1
GENOME="03_genome/genome.fasta"
RAW_DATA_DIR="04_raw_data"
FILT_DIR="05_filtered"
CCS_BAM="$RAW_DATA_DIR/"$SAMPLE".ccs.bam" 
SUBREADS_BAM="$RAW_DATA_DIR/"$SAMPLE".subreads.bam" 

ALIGNED_DIR="06_aligned"

CPU=6

# 0. Create output dir
if [[ ! -d "$ALIGNED_DIR"/pbmm2 ]]
then
  mkdir $ALIGNED_DIR/pbmm2
fi

# 0. Index reference
#pbmm2 index $GENOME 03_genome/genome.fasta.mmi

# 1. Run pbmm2, sort and indem bam output
pbmm2 align $GENOME $CCS_BAM $ALIGNED_DIR/pbmm2/"$SAMPLE".ccs.bam --sort --sample $SAMPLE --bam-index BAI --preset CCS -j $(($CPU / 2)) -J $(($CPU / 2)) --log-level TRACE

pbmm2 align $GENOME $SUBREADS_BAM $ALIGNED_DIR/pbmm2/"$SAMPLE".subreads.bam --sort --sample $SAMPLE --bam-index BAI --preset SUBREAD -j $(($CPU / 2)) -J $(($CPU / 2)) --log-level TRACE

