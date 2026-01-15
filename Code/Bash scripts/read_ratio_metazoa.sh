#!/usr/bin/env bash
# Usage: ./script.sh input_file readcounts_file output_dir

#list of metazoan genomes
infile=$1
#list of readcounts for each metazoan species
readcounts=$2
outfile=$3

#calculates number of reads for each genome for metazoan species

while read -r line; do
#extract genome type
type=$(basename "$line" | awk -F'[_.]' '{print $3}')
#extract species
species=$(basename "$line" | awk -F'[_.]' '{print $1"_"$2}')
#extract total read count for species
readcount=$(awk -v sp="$species" '$2 ~ sp {print $4}' "$readcounts")
#99% chromosome reads
 if [ "$type" = "Chromosome" ]; then
readratio=$(awk -v rc="$readcount" 'BEGIN {print rc * 0.99}')
#1% mitochondria reads
 elif [ "$type" = "Mitochondrion" ]; then
readratio=$(awk -v rc="$readcount" 'BEGIN {print rc * 0.01}')
fi
#append read ratios to the output file
echo -e "$line\t$readratio" >> "$outfile"

done < "$infile"
