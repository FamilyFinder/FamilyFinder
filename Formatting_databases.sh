#!/bin/bash

SCRIPT_DIR=$(dirname "$(readlink -f "$0")")

# Help function
usage() {
  echo "Usage:"
  echo "  formatting_database.sh -i <input> -o <output> [options]"
  echo ""
  echo "Mandatory Arguments:"
  echo "  -i <input>         Input file in FASTA format. The first line must start with '>'."
  echo "  -o <output>        Output file."
  echo ""
  echo "Optional Arguments:"
  echo "  -m <min_length>    Minimum protein length to consider in the database."
  echo "  -M <max_length>    Maximum protein length to consider in the database."
  echo "  -I <domain>        Filter by domain. Specify domains using InterPro classification,"
  echo "                     separated by commas (e.g., IPR001,IPR002,IPR004)."
  echo "  -h                 Display this help message and exit."
  echo ""
  echo "Example:"
  echo "  Formatting_databases.sh -i input.fasta -o output.fasta -m 200 -M 1000 -D IPR001,IPR002,IPR004"
}

while getopts "hi:o:m:M:I:" opt; do
  case ${opt} in
    h )
      usage
      exit 1
      ;;
    i )
      input=$OPTARG
      ;;
    o )
      output=$OPTARG
      ;;
    m )
      min_length=$OPTARG
      ;;
    M )
      max_length=$OPTARG
      ;;
    I )
      domain=$OPTARG
      ;;
    \? )
      echo "Invalid option: $OPTARG" 1>&2
      usage
      exit 1
      ;;
  esac
done

# Check for mandatory arguments
if [ -z "$input" ] || [ -z "$output" ]; then
  echo "Error: Both input and output arguments are mandatory."
  usage
  exit 1
fi


#### Formatting the orthoDB database header
sed -E 's/^>.*"organism_name"[:,]/>/g;s/"pub_gene_id"[:]//g;s/"//g;s/,/_/g;s/description[:]//g;s/ /_/g;s/(}|-)//g;s/[\\]//g;s/[.]//g;s/[;]//g' "$input" > orthodb_header.fasta

# Filter by domain. Important add interproscan in the .bashrc path 
if [[ ! -z "$domain" ]]; 
then
  echo -e "\nRunning interproscan on $domain domain ...\n" 
  interproscan.sh -i orthodb_header.fasta -b interpro_results_database -f tsv 
  [ -d $PWD/temp ] && rm -r $PWD/temp/
  #Python script to keep only the protein with the specific domain in th -D argument
  python3 $SCRIPT_DIR/find_protein_interpro.py interpro_results_database.tsv proteins_to_filter.txt "$domain" 
  # Take only the filtered protein by name
  seqkit grep -f proteins_to_filter.txt orthodb_header.fasta -o filtered_proteins.fasta
  [ -f proteins_to_filter.txt ] && rm proteins_to_filter.txt
  [ -f interpro_results_database.tsv ] && rm interpro_results_database.tsv
else
  mv orthodb_header.fasta filtered_proteins.fasta
fi


# Remove duplicate sequences from the filtered fasta database
seqkit rmdup filtered_proteins.fasta -n -o filtered_proteins_rmdup.fasta

if [ ! -z "$min_length" ] || [ ! -z "$max_length" ]; then
    echo -e "\n Filtering proteins with ${min_length:+Min_length $min_length} ${max_length:+Max_length $max_length} ...\n" 
    seqkit seq ${min_length:+-m $min_length} ${max_length:+-M $max_length} filtered_proteins_rmdup.fasta > "$output"
else
    mv filtered_proteins_rmdup.fasta "$output"
fi

# Remove useless files 
[ -f filtered_proteins.fasta ] && rm filtered_proteins.fasta