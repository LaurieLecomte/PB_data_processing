#!/bin/bash 

# Compare coverage accross mapped bam files

# srun -p medium -c 4 --mem=50G --time=7-00:00:00 -J mosdepth_cov -o log/mosdepth_cov_by_win_%j.log /bin/sh 01_scripts/mosdepth_cov_by_win.sh "safoPUVx_001-21" &

# VARIABLES
SAMPLE=$1
GENOME="03_genome/genome.fasta"
RAW_DATA_DIR="04_raw_data"
FILT_DIR="05_filtered"
ALIGNED_DIR="06_aligned"

SUBREADS_PBMM="$ALIGNED_DIR/pbmm2/"$SAMPLE".subreads.bam"
CCS_PBMM="$ALIGNED_DIR/pbmm2/"$SAMPLE".ccs.bam"

SUBREADS_TRIM="06_aligned/pbmm2/safoPUVx_001-21.subreads.trimmed.bam"

SUBREADS_WIN="$ALIGNED_DIR/winnowmap/"$SAMPLE".subreads.bam"
CCS_WIN="$ALIGNED_DIR/winnowmap/"$SAMPLE".ccs.bam"

COV_DIR="coverage"

WINDOW=1000000
CPU=4

#mosdepth -n --by $WINDOW -t $CPU --thresholds 5,10,20 $COV_DIR/"$SAMPLE".ccs.pbmm2.cov.txt $CCS_PBMM 
mosdepth -n -b $WINDOW -t $CPU $COV_DIR/pbmm2/"$SAMPLE".ccs.pbmm2.cov_win1Mb $CCS_PBMM 

#mosdepth -b $WINDOW -t $CPU --thresholds 1,5,10,20 $COV_DIR/"$SAMPLE".ccs.winnow.cov.txt $CCS_WIN
mosdepth -n -b $WINDOW -t $CPU $COV_DIR/winnowmap/"$SAMPLE".ccs.winnow.cov_win1Mb $CCS_WIN

#mosdepth -n --by $WINDOW -t $CPU --thresholds 1,5,10,20 $COV_DIR/"$SAMPLE".subreads.pbmm2.cov.txt $SUBREADS_PBMM 
mosdepth -n -b $WINDOW -t $CPU $COV_DIR/pbmm2/"$SAMPLE".subreads.pbmm2.cov_win1Mb $SUBREADS_PBMM 

#mosdepth -b $WINDOW -t $CPU --thresholds 1,5,10,20 $COV_DIR/"$SAMPLE".subreads.winnow.cov.txt $SUBREADS_WIN
mosdepth -n -b $WINDOW -t $CPU $COV_DIR/winnowmap/"$SAMPLE".subreads.winnow.cov_win1Mb $SUBREADS_WIN

#mosdepth -n --by $WINDOW -t $CPU --thresholds 1,5,10,20 $COV_DIR/"$SAMPLE".subreads.pbmm2.cov.txt $SUBREADS_TRIM
mosdepth -n -b $WINDOW -t $CPU $COV_DIR/pbmm2/"$SAMPLE".subreads.pbmm2.cov.txt $SUBREADS_TRIM

 