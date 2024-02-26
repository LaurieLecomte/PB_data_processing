#!/bin/bash 

# Produce summary plots for .ccs and subreads files from PacBio Sequel data using longQC
# This script does not work with the longQC version from bioconda, so I needed to install longQC from the GitHub repo and create a custom conda env where I installed python dependencies.
# To use this script:
## 1. Create a custom conda env from longQC_1.2.1_conda_env.txt: conda create --file longQC_1.2.1_conda_env.txt (this contains python dependencies with correct versions)
## 2. Install longQC from the GitHub repo following instructions at https://github.com/yfukasawa/LongQC?tab=readme-ov-file#2-minimap2, in the current directory or somewhere else
## 3. Add the path to the longQC.py script to the LONGQC_EXEC variable below in the script

# To run on all samples in parallel:
# parallel -a 02_infos/ind_PB.txt -j 4 srun -p medium -c 4 --mem=80G --time=2-00:00:00 -J longQC_{} -o log/longQC_{}_%j.log /bin/sh 01_scripts/longQC.sh {} &

# To run on one sample at the time:
# srun -p medium -c 4 --mem=100G --time=7-00:00:00 -J longQC_"safoSABs_027-21" -o log/longQC_"safoSABs_027-21"_%j.log /bin/sh 01_scripts/longQC.sh "safoSABs_027-21" &

# VARIABLES
SAMPLE=$1
RAW_DATA_DIR="04_raw_data"
CCS_BAM="$RAW_DATA_DIR/"$SAMPLE".ccs.bam" 
SUBREAD_BAM="$RAW_DATA_DIR/"$SAMPLE".subreads.bam" 
OUT_DIR="QC"
CPU=4

#LONGQC_EXEC="/project/lbernatchez/users/lalec31/.conda/pkgs/longqc-1.2.0c-hdfd78af_0/bin/longQC.py"
LONGQC_EXEC="/project/lbernatchez/users/lalec31/softwares/LongQC/longQC.py"

# LOAD REQUIRED MODULES

# Create output directory
if [[ ! -d "$OUT_DIR"/longQC/"$SAMPLE" ]]
then
  echo ""$OUT_DIR"/longQC/"$SAMPLE" does not exist"
  mkdir "$OUT_DIR"/longQC/"$SAMPLE"
else 
  echo ""$OUT_DIR"/longQC/"$SAMPLE" exists"
fi


# 1. Run longQC on CCS bam
#python $LONGQC_EXEC sampleqc -x pb-sequel -s $SAMPLE $CCS_BAM -o QC/longQC/$SAMPLE -p $CPU --trim_output QC/longQC/$SAMPLE/trimmed

# 2. Test longQC on subreads bam
# First subset to quicken process
module load samtools/1.15
samtools index $SUBREAD_BAM -@ $CPU
samtools view $SUBREAD_BAM -b -s 0.05 --threads $CPU > $RAW_DATA_DIR/"$SAMPLE".subreads.subs0.05.bam 

python $LONGQC_EXEC sampleqc -x pb-sequel -s $SAMPLE $RAW_DATA_DIR/"$SAMPLE".subreads.subs0.05.bam -o QC/longQC/$SAMPLE -p $CPU --trim_output QC/longQC/$SAMPLE/"$SAMPLE".subreads.subs0.05.trimmed.fastq -d

#python $LONGQC_EXEC sampleqc -x pb-sequel -s $SAMPLE $SUBREAD_BAM -o QC/longQC/$SAMPLE -p $CPU --trim_output QC/longQC/$SAMPLE/"$SAMPLE".subreads.trimmed.fastq -d
#bgzip QC/longQC/$SAMPLE/"$SAMPLE".subreads.trimmed.fastq
