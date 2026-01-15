#!/usr/bin/env bash
# Usage: ./script.sh input_file output_dir

#for all contigs in current directory with _contigs.fasta file ending
#align them to input reference genome and create consensus genome sequence
#input reference genome
ref=$1
#loop over all contigs in currect directory
for f in *_contigs.fasta; do
    base="${f%_contigs.fasta}"
    #align using bowtie2 and sort bam file
    bowtie2 --threads 5 -f -x "$ref" -U "$f" --no-unal | samtools view -q 20 -bS | samtools sort > "${base}.sorted.bam"
    #use angsd to call consensus genome sequence
    angsd -dofasta 2 -docounts 1 -minmapq 20 -minq 25 -i "${base}.sorted.bam" -out "${base}_consensus"
    zcat "${base}_consensus.fa.gz" > "${base}_consensus.fa"
    #rename fasta header to avoid issues with raxml downstream
    sed -i "1s/^>.*/>$base/" "${base}_consensus.fa"
done

