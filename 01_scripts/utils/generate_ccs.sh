#!/bin/bash 

# This script can be used for generating a .ccs file from .subreads file
# It require smrtlink-tools (tested with version 10.1.0), which can be installed in a conda env

# Sample coclEASm_008-20's CCS file we got from sequencing supplier looks unusually small, so we try to generate another .ccs file to see if small size is due to incomplete conversion

# To run on a single sample
# srun -p medium -c 4 --time=7-00:00 -J generate_ccs_coclEASm_008-20 --mem=20G -o log/generate_ccs_coclEASm_008-20_%j.log /bin/sh 01_scripts/generate_ccs.sh "coclEASm_008-20" &

# VARIABLES
SAMPLE=$1
RAW_DATA_DIR="04_raw_data"
SUBREADS_BAM="$RAW_DATA_DIR/"$SAMPLE".subreads.bam" 
CCS_DIR=ccs
CPU=4

# Generate .ccs
ccs --min-passes 3 --num-threads $CPU $SUBREADS_BAM $CCS_DIR/"$SAMPLE".ccs.bam --log-level WARN
