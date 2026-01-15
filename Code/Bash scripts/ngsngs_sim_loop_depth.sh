#!/usr/bin/env bash

#addition to ngsngs_sim_loop.sh. Calculate read depth for each species in simulated metagenomic dataset
#file with percentage of total reads for each species (from bash script read_ratio_percent.sh)
percentfile=$1
#path to references
referencepath=$2
#path to readlength distribution file (from R script rlens.rmd)
lengthdistfile=$3
#number of simulated reads in total
totalreads=$4
outfile=$5


echo -e "Reference\tDepth_of_coverage" > "${outfile}_depth.txt"
#calculate mean read length from length distribution file
meanlength=$(awk '{p=$2-prev; mean+=p*$1; prev=$2} END{print mean}' "$lengthdistfile")
#loop over each reference
while read -r name percent; do
#calculate number of reads
readcount=$(awk -v t="$totalreads" -v p="$percent" 'BEGIN {print int(t*p + 0.5)}')
#calculate bp length of reference
referencelength=$(bioawk -c fastx '{ print length($seq) }' < "$referencepath"/"$name"*)
#calculate read depth
depth=$(awk -v r="$readcount" -v m="$meanlength" -v l="$referencelength" 'BEGIN {print (r*m/l)}')

echo -e "${name}\t${depth}" >> "${outfile}_depth.txt"
done < "$percentfile"
