#!/bin/bash 

# Summary plots for CCS PacBio data using NanoPlot
# parallel -a 02_infos/ind_PB.txt -j 4 srun -p ibis_small -c 2 -J {}_QC_NanoPlot --mem=10G -o log/QC_NanoPlot_{}_%j.log 01_scripts/QC_NanoPlot.sh {} &

# head -n1 02_infos/ind_PB.txt | srun -p ibis_small -c 2 -J "{}"_QC_NanoPlot --mem=10G -o log/QC_NanoPlot_"{}"_%j.log /bin/sh 01_scripts/QC_NanoPlot.sh "{}" &

# Raw data path from valeria : /mnt/ibis/lbernatchez/00_raw_data/2021-09-20_anne-laure_long_reads_coregone_safo_pacbio/

# VARIABLES
SAMPLE=$1
RAW_DATA_DIR="04_raw_data"
CCS_BAM="$RAW_DATA_DIR/"$SAMPLE".ccs.bam" 
SUBREAD_BAM="$RAW_DATA_DIR/"$SAMPLE".subreads.bam" 
OUT_DIR="QC"
CPU=2

# Create output directory
if [[ ! -d "$OUT_DIR" ]]
then
  echo "$OUT_DIR does not exist"
  mkdir $OUT_DIR
else 
  echo "$OUT_DIR exists"
fi

# 1. Run NanoPlot on CCS file
NanoPlot --ubam $CCS_BAM --plots dot --legacy hex --title $SAMPLE -o $OUT_DIR/NanoPlot --N50 --verbose -p "$SAMPLE" --title "$SAMPLE" -t $CPU