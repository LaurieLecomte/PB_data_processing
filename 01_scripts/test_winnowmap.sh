#!/bin/bash 

# Pre-computing high frequency k-mers (e.g., top 0.02% most frequent) in a reference with meryl
# Launch in a conda env where Winnowmap has been installed, or load required modules/wheels


# valeria
# srun -p ibis_large -c 6 --time=21-00:00 -J test_winnowmap --mem=70G -o log/test_winnowmap_%j.log /bin/sh 01_scripts/test_winnowmap.sh "safoPUVx_001-21" &

# VARIABLES
SAMPLE=$1
GENOME="03_genome/genome.fasta"
RAW_DATA_DIR="04_raw_data"
FILT_DIR="05_filtered"
CCS_BAM="$RAW_DATA_DIR/"$SAMPLE".ccs.bam" 
SUBREADS_BAM="$RAW_DATA_DIR/"$SAMPLE".subreads.bam" 

ALIGNED_DIR="06_aligned"

CPU=8

# LOAD REQUIRED MODULES
# module load winnowmap

# 0. Create output dir
if [[ ! -d "$ALIGNED_DIR"/winnowmap ]]
then
  mkdir $ALIGNED_DIR/winnowmap
fi

# O. Convert to bam
# CCS
#pbindex $CCS_BAM
#bam2fasta -o $RAW_DATA_DIR/"$SAMPLE".ccs $CCS_BAM -p "$SAMPLE"

# SUBREADS
#pbindex $SUBREADS_BAM
#bam2fasta -o $RAW_DATA_DIR/"$SAMPLE".subreads $SUBREADS_BAM -p "$SAMPLE"

# 1. Pre-computing high frequency k-mers (e.g., top 0.02% most frequent) in a reference with meryl
#meryl count k=15 output $ALIGNED_DIR/winnowmap/merylDB $GENOME threads=$CPU
#meryl print greater-than distinct=0.9998 $ALIGNED_DIR/winnowmap/merylDB > $ALIGNED_DIR/winnowmap/repetitive_k15.txt


# 2. Run winnowmap on CCS
winnowmap -t $CPU --MD --cs -W $ALIGNED_DIR/winnowmap/repetitive_k15.txt -a -x map-pb -Y -R "@RG\tID:$SAMPLE\tSM:$SAMPLE" $GENOME $RAW_DATA_DIR/"$SAMPLE".ccs.fasta.gz > $ALIGNED_DIR/winnowmap/"$SAMPLE".ccs.sam

# Convert to bam and sort 
samtools view -b $ALIGNED_DIR/winnowmap/"$SAMPLE".ccs.sam --threads $CPU | samtools sort --threads $CPU > $ALIGNED_DIR/winnowmap/"$SAMPLE".ccs.bam
#rm $ALIGNED_DIR/winnowmap/"$SAMPLE".ccs.sam


# 3. Run winnowmap on SUBREADS
winnowmap -t $CPU --MD --cs -W $ALIGNED_DIR/winnowmap/repetitive_k15.txt -a -x map-pb -Y -R "@RG\tID:$SAMPLE\tSM:$SAMPLE" $GENOME $RAW_DATA_DIR/"$SAMPLE".subreads.fasta.gz > $ALIGNED_DIR/winnowmap/"$SAMPLE".subreads.sam

# Convert to bam and sort 
samtools view -b $ALIGNED_DIR/winnowmap/"$SAMPLE".subreads.sam --threads $CPU | samtools sort --threads $CPU > $ALIGNED_DIR/winnowmap/"$SAMPLE".subreads.bam
#rm $ALIGNED_DIR/winnowmap/"$SAMPLE".subreads.sam


# 3. Run winnowmap on TRIMMED SUBREADS
winnowmap -t $CPU --MD --cs -W $ALIGNED_DIR/winnowmap/repetitive_k15.txt -a -x map-pb -Y -R "@RG\tID:$SAMPLE\tSM:$SAMPLE" $GENOME "QC/longQC/$SAMPLE/"$SAMPLE".subreads.trimmed.fastq" > $ALIGNED_DIR/winnowmap/"$SAMPLE".subreads.trimmed.sam

# Convert to bam and sort 
samtools view -b $ALIGNED_DIR/winnowmap/"$SAMPLE".subreads.trimmed.sam --threads $CPU | samtools sort --threads $CPU > $ALIGNED_DIR/winnowmap/"$SAMPLE".subreads.trimmed.bam 
#rm $ALIGNED_DIR/winnowmap/"$SAMPLE".subreads.trimmed.sam

