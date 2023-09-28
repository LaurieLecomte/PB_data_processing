#!/bin/bash

# srun -p small -J plot_coverage -o log/plot_coverage_%j.log /bin/sh 01_scripts/plot_coverage.sh &

# VARIABLES
SAMPLE="safoCHIm_2359-21"
GENOME="03_genome/genome.fasta"

ALIGNED_DIR="06_aligned"
COV_DIR="coverage"

PICARD="/project/lbernatchez/programs/picard/picard.jar"
EXEC="/project/lbernatchez/users/lalec31/softwares/install/jvarkit/dist/wgscoverageplotter.jar"

# LOAD REQUIRED MODULES
module load samtools/1.12


# 0. Create output dir
if [[ ! -d $COV_DIR/pbmm2 ]]
then
  mkdir $COV_DIR/pbmm2
fi 

if [[ ! -d $COV_DIR/winnowmap ]]
then
  mkdir $COV_DIR/winnowmap
fi      

# 0. Index genome
#java -jar $PICARD CreateSequenceDictionary R=$GENOME O="$GENOME".dict

# 1. Plot coverage with winnowmap bam files
## subreads
#samtools index $ALIGNED_DIR/winnowmap/"$SAMPLE".subreads.bam -b
#java -jar $EXEC --dimension 2000x500 -C -1 -R "$GENOME" -I "^CM.+" $ALIGNED_DIR/winnowmap/"$SAMPLE".subreads.bam --percentile median > $COV_DIR/winnowmap/"$SAMPLE".subreads.cov.svg

## ccs

# 2. Plot coverage with pbmm2 bam files
## subreads
samtools index $ALIGNED_DIR/pbmm2/"$SAMPLE".subreads.bam -b 
java -jar $EXEC --dimension 2000x500 -C -1 -R "$GENOME" -I "^CM.+" $ALIGNED_DIR/pbmm2/"$SAMPLE".subreads.bam  --percentile median  > $COV_DIR/pbmm2/"$SAMPLE".subreads.cov.svg

## ccs
#samtools index $ALIGNED_DIR/pbmm2/"$SAMPLE".ccs.bam -b 
#java -jar $EXEC --dimension 2000x500 -C -1 -R "$GENOME" -I "^CM.+" $ALIGNED_DIR/pbmm2/"$SAMPLE".ccs.bam  --percentile median  > $COV_DIR/pbmm2/"$SAMPLE".ccs.cov


