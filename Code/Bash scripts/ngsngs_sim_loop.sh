#!/usr/bin/env bash

#create simulated metagenomic dataset
#file with percent reads of for each reference (from bash script read_ratio_percent.sh)
percentfile=$1
#path to references for simulation (from bash script get_references.sh)
referencepath=$2
#read length distribution file (from R script rlens.rmd)
lengthdistfile=$3
#nucleotide substitution file (from R script Substitution_rate.rmd)
nucmismatchfile=$4
#number of reads to simulate in total
totalreads=$5
#output files names
outfile=$6

#loop over species names and percentages
while read -r name percent; do
#defines number of reads to create for each reference using percentage file
readcount=$(awk -v t="$totalreads" -v p="$percent" 'BEGIN {print int(t*p + 0.5)}')
#simulate reads using ngsngs
ngsngs -i "$referencepath"/"$name" -r "$readcount" -lf "$lengthdistfile" -f fq -qs 40 -seq PE -mf "$nucmismatchfile" -o "${outfile}_${name}"

#concatenate the simulated reads in a single file
cat "${outfile}_${name}_R1.fq" >> "${outfile}_all_R1.fq"
cat "${outfile}_${name}_R2.fq" >> "${outfile}_all_R2.fq"
done < "$percentfile"
