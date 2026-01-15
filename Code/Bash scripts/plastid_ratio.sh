#!/usr/bin/env bash
# Usage: ./script.sh input_file output_dir

#calculate plastid read ratios for Tjornin lake samples
infile=$1
outfile=$2

#create headers for outfile
echo -e "Sample\tCompetitively_Mapped_Viridiplantae_Reads\tChloroplast_Mapped_Reads\tChloroplast_to_Viridiplantae_ratio" > "$outfile"

#for each sample, check for both chloroplast alignment and competitive alignment
while IFS= read -r statfile; do
#take sample names from statfile of competitive alignments
 sample=$(basename "$statfile" | awk -F'_' '{print $1"_"$2}')
 #check for sample names in directory of chloroplast alignments
 bamfile="/projects/wintherpedersen/scratch/august_chloroplast/"$sample"_collapsed_full.bam"
 if [ ! -f "$bamfile" ]; then
        echo "Warning: BAM file not found for $sample, skipping."
        continue
    fi

#get number of reads for Viridiplantae from stats file from competitive mapping
 CompPlantReads=$(zcat "$statfile" | awk '$2 ~ /Viridiplantae/ {print $5}')
#get mapped reads from chloroplast aligned bam file
 PlastidReads=$(samtools stats "$bamfile" | awk '/^SN/ && /reads mapped:/{print $NF}')
#calculate ratio of chloroplast reads to competively mapped viridiplantae reads
 PlastidToTotal=$(echo "scale=6; $PlastidReads / $CompPlantReads" | bc)
#append values to output file
 echo -e "$sample\t$CompPlantReads\t$PlastidReads\t$PlastidToTotal" >> "$outfile"

done < "$infile"
