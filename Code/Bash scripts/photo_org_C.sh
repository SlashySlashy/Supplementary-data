#!/usr/bin/env bash
# Usage: ./script.sh input_directory output_dir

#make file with read counts of 4 taxa extracted from files in directory
indirectory=$1
outfile=$2

echo -e "CGG_ID\tViridiplantae\tRhodophyta\tOchrophyta\tCyanobacteriota" > "$outfile"

#input directory must contain files where 1st column is taxid and 5th is number of reads
for file in "$indirectory"/*; do
 if [ -f "$file" ]; then
 #extract sample name
 CGG=$(basename "$file" | awk -F'.' '{print $1}')
 #viridiplantae reads
 Vreads=$(zcat "$file" | awk '$1==33090 {print $5}')
 #rhodophyta reads
 Rreads=$(zcat "$file" | awk '$1==2763 {print $5}')
 #ochrophyta reads
 Oreads=$(zcat "$file" | awk '$1==2696291 {print $5}')
 #cyanobacteriota reads
 Creads=$(zcat "$file" | awk '$1==1117 {print $5}')
 #add sample name and read counts to output file
 echo -e "$CGG\t$Vreads\t$Rreads\t$Oreads\t$Creads" >> "$outfile"
 fi
done
