#!/bin/bash 

# Align PacBio reads with winnowmap
# Launch in a conda env where winnowmap 2.02, smrtlink-tools 10.1.0 and samtools 1.15 have been installed, or load required modules/wheels

# manitou
# parallel -a 02_infos/ind_PB.txt -j 4 srun -p medium -c 4 --time=7-00:00 -J 02_winnowmap_{} --mem=70G -o log/02_winnowmap_{}_%j.log /bin/sh 01_scripts/02_winnowmap.sh {} &

# valeria
# parallel -a 02_infos/ind_PB.txt -j 4 srun -p ibis_medium -c 6 --time=7-00:00 -J 02_winnowmap_{} --mem=70G -o log/02_winnowmap_{}_%j.log /bin/sh 01_scripts/02_winnowmap.sh {} &

# VARIABLES
SAMPLE=$1
GENOME="03_genome/genome.fasta"
RAW_DATA_DIR="04_raw_data"
FILT_DIR="05_filtered"
CCS_BAM="$RAW_DATA_DIR/"$SAMPLE".ccs.bam" 
SUBREADS_BAM="$RAW_DATA_DIR/"$SAMPLE".subreads.bam" 

ALIGNED_DIR="06_aligned"

CPU=4

# LOAD REQUIRED MODULES
#module load samtools/1.15

# 0. Create output dir
if [[ ! -d "$ALIGNED_DIR"/winnowmap ]]
then
  mkdir $ALIGNED_DIR/winnowmap
fi

# O. Convert .ccs bam file to fasta
pbindex $CCS_BAM
bam2fasta -o $RAW_DATA_DIR/"$SAMPLE".ccs $CCS_BAM -p "$SAMPLE"


# 1. Run winnowmap on converted file, using pre-computed kmers
## -Y enables mapping of soft-clipped reads 
winnowmap -t $CPU --MD --cs -W $ALIGNED_DIR/winnowmap/repetitive_k15.txt -a -x map-pb -Y -R "@RG\tID:$SAMPLE\tSM:$SAMPLE" $GENOME $RAW_DATA_DIR/"$SAMPLE".ccs.fasta.gz > $ALIGNED_DIR/winnowmap/"$SAMPLE".ccs.sam

# Convert to bam, sort and index 
samtools sort --threads $CPU  $ALIGNED_DIR/winnowmap/"$SAMPLE".ccs.sam -o $ALIGNED_DIR/winnowmap/"$SAMPLE".ccs.bam
samtools index $ALIGNED_DIR/winnowmap/"$SAMPLE".ccs.bam


# Clean up
#rm $ALIGNED_DIR/winnowmap/"$SAMPLE".ccs.sam
