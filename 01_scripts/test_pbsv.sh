#!/bin/bash 

# Test pbsv using mapped reads from pbmm2

# srun -p medium -c 4 --time=7-00:00 -J pbsv_safoCHIm_2359-21 --mem=20G -o log/pbsv_safoCHIm_2359-21_%j.log /bin/sh 01_scripts/test_pbsv.sh "safoCHIm_2359-21" &

# VARIABLES
SAMPLE=$1
GENOME="03_genome/genome.fasta"
RAW_DATA_DIR="04_raw_data"
FILT_DIR="05_filtered"
CCS_BAM="$RAW_DATA_DIR/"$SAMPLE".ccs.bam" 
SUBREADS_BAM="$RAW_DATA_DIR/"$SAMPLE".subreads.bam" 

ALIGNED_DIR="06_aligned"

MAPPED_SUB="$ALIGNED_DIR/pbmm2/"$SAMPLE".subreads.bam"
MAPPED_CCS="$ALIGNED_DIR/pbmm2/"$SAMPLE".ccs.bam"
#MAPPED_TRIM_SUB
#MAPPED_TRIM_CCS

SV_DIR="pbsv_test"

CPU=4

# Discover signatures of structural variation
pbsv discover $MAPPED_SUB $SV_DIR/"$SAMPLE".subreads.svsig.gz --log-level WARN

# optionally index svsig.gz to allow random access via `pbsv call -r`
tabix -c '#' -s 3 -b 4 -e 4 $SV_DIR/"$SAMPLE".subreads.svsig.gz

# It is highly recommended to provide one tandem repeat annotation .bed file of your reference to pbsv discover via --tandem-repeats. This increases sensitivity and recall. Feel free to use the following for human SV calling: GRCh38 or hs37d5/hg19.

#Sample names are transferred from the RG headers to the .svsig.gz file.

# Call structural variants and assign genotypes
# Call structural variants from structural variant signatures, jointly for all samples of interest. One or more .svsig.gz files are accepted, including multiple .svsig.gz for a single sample and/or svsig.gz for multiple samples. If the input is CCS reads, please add --ccs to the following call:
#pbsv call ref.fa ref.sample1.svsig.gz ref.sample2.svsig.gz ref.var.vcf

pbsv call $GENOME $SV_DIR/"$SAMPLE".subreads.svsig.gz $SV_DIR/"$SAMPLE".subreads.vcf -j $CPU --log-level WARN

# CCS
#pbsv discover $MAPPED_CCS $SV_DIR/"$SAMPLE".ccs.svsig.gz --log-level INFO -s $SAMPLE
#tabix -c '#' -s 3 -b 4 -e 4 $SV_DIR/"$SAMPLE".ccs.svsig.gz
#pbsv call $GENOME $SV_DIR/"$SAMPLE".ccs.svsig.gz $SV_DIR/"$SAMPLE".ccs.vcf --ccs -j $CPU --log-level INFO


# From winnowmap
pbsv discover "$ALIGNED_DIR/winnowmap/"$SAMPLE".subreads.bam" $SV_DIR/"$SAMPLE".subreads.win.svsig.gz --log-level INFO -s $SAMPLE
tabix -c '#' -s 3 -b 4 -e 4 $SV_DIR/"$SAMPLE".subreads.win.svsig.gz
pbsv call $GENOME $SV_DIR/"$SAMPLE".subreads.win.svsig.gz $SV_DIR/"$SAMPLE".subreads.win.vcf -j $CPU --log-level INFO


pbsv discover "$ALIGNED_DIR/winnowmap/"$SAMPLE".ccs.bam" $SV_DIR/"$SAMPLE".ccs.win.svsig.gz --log-level INFO -s $SAMPLE
tabix -c '#' -s 3 -b 4 -e 4 $SV_DIR/"$SAMPLE".ccs.win.svsig.gz
pbsv call $GENOME $SV_DIR/"$SAMPLE".ccs.win.svsig.gz $SV_DIR/"$SAMPLE".ccs.win.vcf --ccs -j $CPU --log-level INFO