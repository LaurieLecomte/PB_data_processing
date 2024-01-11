#!/bin/bash 

# Summary plots for PacBio data using longQC
# parallel -a 02_infos/ind_PB.txt -j 4 /project/lbernatchez/users/lalec31/projets_labo/FISHES/02_long_reads/PB_data_processing/PB_data_processing/01_scripts/

# head -n2 02_infos/ind_PB.txt | srun -p ibis_small -c 2 -J "{}"_QC_longQC --mem=10G -o log/QC_longQC_"{}"_%j.log /bin/sh 01_scripts/QC_longQC.sh "{}" &

# Raw data path from valeria : /mnt/ibis/lbernatchez/00_raw_data/2021-09-20_anne-laure_long_reads_coregone_safo_pacbio/

# srun -p medium -c 8 --mem=100G --time=7-00:00:00 -J "safoSABs_027-21"_QC_longQC -o log/QC_longQC_"safoSABs_027-21"_%j.log /bin/sh 01_scripts/QC_longQC.sh "safoSABs_027-21" &

# VARIABLES
SAMPLE=$1
RAW_DATA_DIR="04_raw_data"
CCS_BAM="$RAW_DATA_DIR/"$SAMPLE".ccs.bam" 
SUBREAD_BAM="$RAW_DATA_DIR/"$SAMPLE".subreads.bam" 
OUT_DIR="QC"
CPU=8

#LONGQC_EXEC="/project/lbernatchez/users/lalec31/.conda/pkgs/longqc-1.2.0c-hdfd78af_0/bin/longQC.py"
LONGQC_EXEC="/project/lbernatchez/users/lalec31/softwares/LongQC/longQC.py"

# LOAD REQUIRED MODULES
#module load python/3.9

# Create output directory
#if [[ ! -d "$OUT_DIR"/longQC/"$SAMPLE" ]]
#then
#  echo ""$OUT_DIR"/longQC/"$SAMPLE" does not exist"
#  mkdir "$OUT_DIR"/longQC/"$SAMPLE"
#else 
#  echo ""$OUT_DIR"/longQC/"$SAMPLE" exists"
#fi


# 1. Run longQC on CCS bam
#python $LONGQC_EXEC sampleqc -x pb-sequel -s $SAMPLE $CCS_BAM -o QC/longQC/$SAMPLE -p $CPU --trim_output QC/longQC/$SAMPLE/trimmed

# Run longQC on subreads bam
# First subset to quicken process
module load samtools/1.15
samtools index $SUBREAD_BAM -@ $CPU
samtools view $SUBREAD_BAM -b -s 0.05 --threads $CPU > $RAW_DATA_DIR/"$SAMPLE".subreads.subs0.05.bam 

python $LONGQC_EXEC sampleqc -x pb-sequel -s $SAMPLE $RAW_DATA_DIR/"$SAMPLE".subreads.subs0.05.bam -o QC/longQC/$SAMPLE -p $CPU --trim_output QC/longQC/$SAMPLE/"$SAMPLE".subreads.subs0.05.trimmed.fastq -d

#python $LONGQC_EXEC sampleqc -x pb-sequel -s $SAMPLE $SUBREAD_BAM -o QC/longQC/$SAMPLE -p $CPU --trim_output QC/longQC/$SAMPLE/"$SAMPLE".subreads.trimmed.fastq -d
#bgzip QC/longQC/$SAMPLE/"$SAMPLE".subreads.trimmed.fastq
