#!/usr/bin/env bash
# Usage: ./script.sh input_file readcounts_file output_dir

#file with viridiplantae genomes
infile=$1
#file with read counts for viridiplantae species
readcounts=$2
outfile=$3

#calculate the number of reads for each genome for each species
#loop over viridiplantae genomes
while read -r line; do
#extract genome type
type=$(basename "$line" | awk -F'[_.]' '{print $3}')
#extract species name
species=$(basename "$line" | awk -F'[_.]' '{print $1"_"$2}')
#extract combined read count for the species
readcount=$(awk -v sp="$species" '$2 ~ sp {print $4}' "$readcounts")
#94.555% chromosome reads
 if [ "$type" = "Chromosome" ]; then
readratio=$(awk -v rc="$readcount" 'BEGIN {print rc * 0.94555}')
#1% mitochondria reads
 elif [ "$type" = "Mitochondrion" ]; then
readratio=$(awk -v rc="$readcount" 'BEGIN {print rc * 0.01}')
#4.445% chloroplast reads
 elif [ "$type" = "Chloroplast" ]; then
readratio=$(awk -v rc="$readcount" 'BEGIN {print rc * 0.04445}')
fi
echo -e "$line\t$readratio" >> "$outfile"

done < "$infile"
