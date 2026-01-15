#!/usr/bin/env bash
# Usage: ./script.sh input_file output_dir

#converts output of read_ratio_viridiplantae.sh and read_ratio_metazoa.sh from read counts to percentages of total
#infile is concatenated output of read_ratio_viridiplantae.sh and read_ratio_metazoa.sh
infile=$1
outfile=$2

#calculate total sum of read counts
sum=$(awk '{s+=$2}END{print s}' "$infile")
#loop over each genome
while read -r name value; do
#calculate percentage of total read count for the species
percent=$(awk -v s="$sum" -v v="$value" 'BEGIN {print v/s}')
echo -e "$name\t$percent" >> "$outfile"
done < "$infile"
