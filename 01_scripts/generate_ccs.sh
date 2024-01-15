#!/bin/bash 

# Sample coclEASm_008-20's CCS file we got from sequencing supplier looks unusually small, so we try to generate another CCS file

# srun -p medium -c 4 --time=7-00:00 -J generate_ccs_coclEASm_008-20 --mem=20G -o log/generate_ccs_coclEASm_008-20_%j.log /bin/sh 01_scripts/generate_ccs.sh "coclEASm_008-20" &

# VARIABLES
SAMPLE=$1
RAW_DATA_DIR="04_raw_data"
SUBREADS_BAM="$RAW_DATA_DIR/"$SAMPLE".subreads.bam" 

CCS_DIR=ccs

CPU=4


#CL /cvmfs/soft.mugqic/CentOS6/software/SMRTLink/SMRTLink-11.0.0/install/smrtlink-release_11.0.0.146107/bundles/smrttools/install/smrttools-release_11.0.0.146107/private/pacbio/unanimity/binwrap/../../../../private/pacbio/unanimity/bin/ccs --min-passes 3 --num-threads 20 /lustre03/project/6033481/4nanuq2/platform/sequelRuns/r64128_20220825_181935/1_A01/m64128_220825_183034.subreads.bam /lustre03/project/6033481/4nanuq2/platform/sequelRuns/r64128_20220825_181935/1_A01/m64128_220825_183034.ccs.bam

ccs --min-passes 3 --num-threads $CPU $SUBREADS_BAM $CCS_DIR/"$SAMPLE".ccs.bam --log-level WARN
