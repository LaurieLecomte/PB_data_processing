#!/bin/bash 

# Produce summary plots for .ccs files from PacBio Sequel data using longQC. Because we use .ccs files, adaptors should already have been trimmed, so we won't save the trimmed read output. LongQc also does not understand ccs quality scores, so quality metrics based on read quality are not useful here. 

# This script does not work with the longQC version from bioconda, so I needed to install longQC from the GitHub repo and create a custom conda env where I installed python dependencies.
# To use this script:
## 1. Create a custom conda env from longQC_1.2.1_conda_env.txt: conda create --file longQC_1.2.1_conda_env.txt (this contains python dependencies with correct versions)
## 2. Install longQC from the GitHub repo following instructions at https://github.com/yfukasawa/LongQC?tab=readme-ov-file#2-minimap2, in the current directory or somewhere else
## 3. Add the path to the longQC.py script (located in the installation directory) to the LONGQC_EXEC variable below in the script, or save and export as an env variable in the current session (LONGQC_EXEC="/project/lbernatchez/users/lalec31/softwares/LongQC/longQC.py" then export LONGQC_EXEC)

# To run on all samples in parallel:
# parallel -a 02_infos/ind_PB.txt -j 4 srun -p medium -c 6 --mem=80G --time=2-00:00:00 -J longQC_ccs_{} -o log/longQC_ccs_{}_%j.log /bin/sh 01_scripts/longQC_ccs.sh {} &

# To run on one sample at the time:
# srun -p medium -c 6 --mem=80G --time=2-00:00:00 -J longQC_ccs_"safoSABs_027-21" -o log/longQC_ccs_"safoSABs_027-21"_%j.log /bin/sh 01_scripts/longQC_ccs.sh "safoSABs_027-21" &

# VARIABLES
SAMPLE=$1
RAW_DATA_DIR="04_raw_data"
CCS_BAM="$RAW_DATA_DIR/"$SAMPLE".ccs.bam" 
SUBREAD_BAM="$RAW_DATA_DIR/"$SAMPLE".subreads.bam" 
OUT_DIR="QC"
CPU=6

#LONGQC_EXEC="/project/lbernatchez/users/lalec31/.conda/pkgs/longqc-1.2.0c-hdfd78af_0/bin/longQC.py"
#LONGQC_EXEC="/project/lbernatchez/users/lalec31/softwares/LongQC/longQC.py"

# LOAD REQUIRED MODULES

# Create output directory
if [[ ! -d "$OUT_DIR"/longQC ]]
then
  echo ""$OUT_DIR"/longQC does not exist"
  mkdir "$OUT_DIR"/longQC
else 
  echo ""$OUT_DIR"/longQC exists"
fi



# 1. Run longQC on .ccs bam
python $LONGQC_EXEC sampleqc -x pb-sequel -s $SAMPLE $CCS_BAM -o QC/longQC/$SAMPLE -p $CPU -d


