#!/bin/bash

set -e  # Exit if any command fails
set -x	# Print commands as they execute

# Input reference fasta, read lengths distribution, and read depth of coverage
REFERENCE=$1
COVERAGE=$3
LENGTH=$2

#Names for output files
BASE_NAME="${REFERENCE%.*}_${COVERAGE}_depth"

#Use NGS to simulate reads
NGSNGS/ngsngs -i ./"$REFERENCE" -c "$COVERAGE" -lf "$LENGTH" -f fq -qs 40 -seq PE -circ -mf ./nucleotide_mismatch_table.txt -o "$BASE_NAME"

#MEGAhit assembly
megahit -1 ${BASE_NAME}_R1.fq -2 ${BASE_NAME}_R2.fq -o Megahit_$BASE_NAME

#Rename and move megahit contigs
mv Megahit_${BASE_NAME}/final.contigs.fa Megahit_"$BASE_NAME".fa

#Merge reads for CarpeDeam assembly
flash ${BASE_NAME}_R1.fq ${BASE_NAME}_R2.fq -o ${BASE_NAME}_merged

#CarpeDeam assembly
carpedeam ancient_assemble ${BASE_NAME}_merged.extendedFrags.fastq CarpeDeam_$BASE_NAME.fasta tmp --ancient-damage d_high

#SPAdes assembly
SPAdes-4.0.0-Linux/bin/spades.py -1 ${BASE_NAME}_R1.fq -2 ${BASE_NAME}_R2.fq --phred-offset 33 -o SPAdes_$BASE_NAME

#Rename and move SPAdes contigs
mv SPAdes_${BASE_NAME}/contigs.fasta SPAdes_"$BASE_NAME".fasta

echo "Workflow completed successfully."
