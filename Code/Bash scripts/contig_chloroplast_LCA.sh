#!/usr/bin/env bash
# Usage: ./script.sh contig_file ref_directory output_name

#create LCA from contigs aligned against chloroplast references
CONTIG=$1
REF=$2
OUTPUT=$3

#uses minimap2 to align contigs against reference genomes
minimap2 -a "$REF" "$CONTIG" -o "${OUTPUT}_minimap2.sam"
#converts sam file to bam and sort it
samtools view -S -b "${OUTPUT}_minimap2.sam" | samtools sort -n > "${OUTPUT}_minimap2.sorted.bam"
rm "${OUTPUT}_minimap2.sam"
#for each contig, removes all alignments with more mismatches than the best alignment
python best_as_filter.py "${OUTPUT}_minimap2.sorted.bam" "${OUTPUT}_minimap2.sorted.filtered.bam"

#use ngsLCA to find LCA taxid for each contig
ngsLCA -names /datasets/caeg_dataset/references/ncbi/20250530/taxonomy/ncbi/taxdump/names.dmp -nodes /datasets/caeg_dataset/references/ncbi/20250530/taxonomy/ncbi/taxdump/nodes.dmp -acc2tax <(zcat /datasets/caeg_dataset/references/ncbi/20250530/taxonomy/ncbi/*.acc2taxid.gz) -simscorelow 0.90 -simscorehigh 1 -fix_ncbi 0 -bam "${OUTPUT}_minimap2.sorted.filtered.bam" -outnames "${OUTPUT}_minimap2.sorted.filtered"

#finds number of assigned contigs for each taxid
cut -f 4 -d : "${OUTPUT}_minimap2.sorted.filtered.lca" | awk '{print $2}' | tail -n +3 | sort | uniq -c | sort -k1 -rn > "${OUTPUT}_minimap2.taxidcounts"
#finds total length of assigned contigs for each taxid
cut -f 3,4,6 -d : "${OUTPUT}_minimap2.sorted.filtered.lca" | awk -F':' '{print $1, $2, $3}' | tail -n +3 | awk '{print $1, $3, $4}' > "${OUTPUT}_minimap2.lengths"
