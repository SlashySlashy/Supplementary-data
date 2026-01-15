#!/usr/bin/env bash
# Usage: ./script.sh input_file output_dir

#script for downloading reference genomes from ncbi from list of species names, using datasets command line tool
#input file with list of taxids and species name and output directory
infile=$1
outdir=$2

mkdir -p "$outdir"


while IFS=$'\t' read -r taxa_ID taxa_name _ _; do
    [[ "$taxa_ID" == "taxa_ID" ]] && continue
    echo "Processing: $taxa_name"

# Check for any reference genomes
if datasets summary genome taxon "$taxa_ID" --reference --as-json-lines \
  | jq -s -e 'length > 0' >/dev/null; then
  echo "Reference genome found. Downloading sequences..."
else
  echo "No reference genome found for $taxa_name" >&2
  continue
fi

# Download all reference genomes for taxid
datasets download genome taxon "$taxa_ID" \
  --reference \
  --include genome,seq-report \
  --filename "${taxa_ID}_ref.zip"


    unzip -p "${taxa_ID}_ref.zip" "ncbi_dataset/data/*/sequence_report.jsonl" > "${taxa_ID}_sequence_report.jsonl"
    fna=$(unzip -Z1 "${taxa_ID}_ref.zip" | grep '_genomic\.fna$')

#loops over mitchondria, chloroplast, and chromosomes to select only one of each as reference, to avoid duplicate references
for type in Mitochondrion Chloroplast Chromosome; do
    ids_file="${taxa_ID}_${type}_acc.txt"
    fasta_out="${outdir}/${taxa_name}_${type}.fasta"

#select the first complete reference for each type and save header to id file
    awk -v type="$type" '
      $0 ~ "\"role\":\"assembled-molecule\"" && $0 ~ "\"assignedMoleculeLocationType\":\""type"\"" {
        if (match($0, /"genbankAccession":"([^"]+)"/, m)) {
          print m[1]
        }
	if (match($0, /"refseqAccession":"([^"]+)"/, n)) {
          print n[1]
        }
      }
    ' "${taxa_ID}_sequence_report.jsonl" > "$ids_file"

#extract reference headers saved in id file from complete reference fasta, save to new fasta file for each organelle
    if [[ -s $ids_file ]]; then
      : > "$fasta_out"
      unzip -p "${taxa_ID}_ref.zip" "$fna" | awk -v list="$ids_file" '
  BEGIN {
    while ((getline line < list) > 0) {
      sub(/\r$/,"",line); if (line!="") wanted[line]=1
    }
  }
   /^>/ {
    print_flag = 0
    for (id in wanted) {
      if (index($0, id) > 0) { print_flag = 1; break }
    }
  }
  { if (print_flag) print } ' > "$fasta_out"

      if [[ -s "$fasta_out" ]]; then
        echo "Saved $fasta_out"
      else
        echo "No ${type} sequences matched headers in FASTA for $taxa_name"
        rm -f "$fasta_out"
      fi
    else
      echo "No ${type} accessions listed for $taxa_name"
    fi
  done

  # Cleanup
  rm -f "${taxa_ID}_ref.zip" "${taxa_ID}_sequence_report.jsonl" "${taxa_ID}_fna.list" "${taxa_ID}_"*_acc.txt


done < <(tr -d '\r' < "$infile")
