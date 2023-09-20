#!/bin/bash 

# Summary plots for CCS PacBio data using longQC
# parallel -a 02_infos/ind_PB.txt -j 4 srun -p ibis_small -c 2 -J {}_QC_longQC --mem=10G -o log/QC_longQC_{}_%j.log 01_scripts/QC_longQC.sh {} &

# head -n2 02_infos/ind_PB.txt | srun -p ibis_small -c 2 -J "{}"_QC_longQC --mem=10G -o log/QC_longQC_"{}"_%j.log /bin/sh 01_scripts/QC_longQC.sh "{}" &

# Raw data path from valeria : /mnt/ibis/lbernatchez/00_raw_data/2021-09-20_anne-laure_long_reads_coregone_safo_pacbio/

# VARIABLES
SAMPLE=$1
RAW_DATA_DIR="04_raw_data"
CCS_BAM="$RAW_DATA_DIR/"$SAMPLE".ccs.bam" 
SUBREAD_BAM="$RAW_DATA_DIR/"$SAMPLE".subreads.bam" 
OUT_DIR="QC"
CPU=2

LONGQC_EXEC="/project/lbernatchez/users/lalec31/.conda/pkgs/longqc-1.2.0c-hdfd78af_0/bin/longQC.py"

# Create output directory
if [[ ! -d "$OUT_DIR"/longQC/"$SAMPLE" ]]
then
  echo ""$OUT_DIR"/longQC/"$SAMPLE" does not exist"
  mkdir "$OUT_DIR"/longQC/"$SAMPLE"
else 
  echo ""$OUT_DIR"/longQC/"$SAMPLE" exists"
fi


# 1. Run longQC on CCS bam
python $LONGQC_EXEC sampleqc -x pb-sequel -s $SAMPLE $CCS_BAM -o QC/longQC/$SAMPLE -p $CPU --trim_output QC/longQC/$SAMPLE/trimmed


