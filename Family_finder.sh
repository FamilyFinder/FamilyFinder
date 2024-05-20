#!/bin/bash

#SBATCH -n 20
#SBATCH -p long
#SBATCH --mem=80000
#SBATCH -t 05-20:00:00
#SBATCH -J anotación
#SBATCH -e /mnt/lustre/home/fblasio/error.err
#SBATCH -o /mnt/lustre/home/fblasio/output

SCRIPT_DIR=$(dirname "$(readlink -f "$0")")

# Help function
usage() {
  echo "Usage:"
  echo "  Family_finder.sh -i <input> -o <output> -d <database> [options]"
  echo ""
  echo "Mandatory Arguments:"
  echo "  -i <input>         Path to the directory containing the genomes."
  echo "  -d <database>      Database for protein predictions (e.g., orthoDB). This database"
  echo "                     will be used with miniprot (https://github.com/lh3/miniprot)."
  echo ""
  echo "Optional Arguments:"
  echo "  -o <output>        Output directory (default: ./Family_finder_out)."
  echo "  -D <database>      Secondary database for predictions, not used with miniprot."
  echo "                     Used alongside -d database by metaeuk."
  echo "  -m <min_length>    Minimum protein length to consider after predictions."
  echo "  -M <max_length>    Maximum protein length to consider after predictions."
  echo "  -H <file.hmm>      HMM file to filter prediction results (recommended when using -D)."
  echo "  -t                 Specify to generate a phylogenetic tree after predictions with IQtree."
  echo "  -q <iqtree_params> Specify iqtree parameters (default: -bb 1000 -alrt 1000 -abayes -mem 20G -safe -T AUTO)."
  echo "  -s                 Specify to retrieve species names."
  echo "                     Note: This may significantly increase processing time."
  echo "  -I <domain>        Specify domains using InterPro classification,"
  echo "                     separated by commas (e.g., IPR001,IPR002,IPR004)."
  echo "  -e <environment>   Specify a Conda environment for third-party components."
  echo "  -h                 Display this help message and exit."
  echo ""
  echo "Example:"
  echo "  Family_finder.sh -i /path/to/genomes -d path/to/orthoDB -o /path/to/output"
}

# defaults
input=""
output="$(realpath ./Family_finder_out)"
min_length=""
max_length=""
perfil_HMM=""
iqtree=false
species_name=false
database=""
iq_params=" -bb 1000 -alrt 1000 -abayes -mem 20G -safe -T AUTO"
database_secondary=""
domain=""
env=""
GREEN='\033[0;32m'
NC='\033[0m' # No Color
RED='\033[31m'

while getopts ":i:o:m:M:H:tq:sd:D:I:e:h" opt; do
  case ${opt} in
    i)
      input="$(realpath $OPTARG)"
      ;;
    o)
      output="$(realpath $OPTARG)"
      ;;
    m)
      min_length="$OPTARG"
      ;;
    M)
      max_length="$OPTARG"
      ;;
    H)
      perfil_HMM="$OPTARG"
      ;;
    t)
      iqtree=true
      ;;
    q)
      iq_params="$OPTARG"
      ;;
    s)
      species_name=true
      ;;
    d)
      database="$(realpath $OPTARG)"
      ;;
    D)
      database_secondary="$(realpath $OPTARG)"
      ;;
    I)
      domain="$OPTARG"
      ;;
    e)
      env="$OPTARG"
      ;;
    h)
      usage
      exit 0
      ;;
    \?)
      echo "Invalid option: -$OPTARG. The following arguments are required: -i (input) and -d (database)." >&2
      usage
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

shift $((OPTIND-1))

################# CHECK IF THE INPUT PARAMETERS ARE OK #################

# Check for mandatory parameters
if [[ -z "$input" || -z "$database" ]]; then
  echo "Missing required arguments: -i (input) and -d (database) are mandatory." >&2
  usage
  exit 1
fi

# Verify input files and directories
if [[ ! -d "$input" ]]; then
  echo "Input directory does not exist: $input" >&2
  exit 1
fi

if [[ ! -f "$database" ]]; then
  echo "Database file does not exist: $database" >&2
  exit 1
fi

if [[ ! -z "$database_secondary" && ! -f "$database_secondary" ]]; then
  echo "Secondary database file does not exist: $database_secondary" >&2
  exit 1
fi

# Check if the database files are in FASTA format
if ! grep -q "^>" "$database"; then
  echo "The database file is not in FASTA format: $database" >&2
  exit 1
fi

# Check if the secondary database is in FASTA format
if [[ ! -z "$database_secondary" ]]; then
    if ! grep -q "^>" "$database_secondary"; then
        echo "Error: Check your $database_secondary is a valid FASTA file" >&2
        exit 1
    fi
fi

# Check if the directory exist or not. If the output directory doesn't exist make one
if [ ! -d "$output" ]; then
    mkdir -p "$output"
fi

echo -e "###### LOG ######\n\n" > "$output"/program.log

# Activate Conda environment if specified
if [[ ! -z "$env" ]]; then
  # Source the Conda base script
  eval "$(conda shell.bash hook)"
  CONDA_BASE=$(conda info --base);
  source "${CONDA_BASE}/etc/profile.d/conda.sh";

  # Activate the environment
  conda activate "$env";
  echo -e "\nConda environment $env ${GREEN}activated${NC}\n" 
  echo -e "\nConda environment $env activated\n" >> "$output"/program.log
fi

########################### START ###########################
# Store variable 
input_folder=($(find "$input" -type f \( -name "*.fna" -o -name "*.fasta" \)))

database_file="$database"



echo -e "Checking dependencies ..."

# CHECK dependencies
echo -e "Checking mandatory dependencies ... \n" >> "$output"/program.log

dependencies=(seqkit python3 miniprot Rscript metaeuk bedtools)

if [[ "$species_name" = true || ! -z "$domain" ]]; then
    dependencies+=(interproscan.sh esearch)
fi

if [[ ! -z "$perfil_HMM" ]]; then
    dependencies+=(hmmbuild)
fi

if [[ "$iqtree" = true ]]; then
    dependencies+=(iqtree kalign trimal)
fi

for program in "${dependencies[@]}"
do
  if ! command -v $program &>/dev/null; then
    echo -e "\n$program needs to be installed or check that the program is in your PATH variable." >> "$output"/program.log
    echo -e "\n$program needs to be installed or check that the program is in your PATH variable.\n"
    exit 1
  else
    echo -e "\n$program ${GREEN}OK${NC}" 
    echo -e "\n$program OK" >> "$output"/program.log
  fi
done

echo -e "\nChecking R dependencies ...\n"
# Check R packages
if ! Rscript -e 'packages <- c("dplyr", "tidyr", "writexl", "readxl", "stringr", "data.table"); for (pkg in packages) {if(!suppressWarnings(suppressPackageStartupMessages(suppressMessages(require(pkg, character.only = TRUE))))) {cat(paste0(pkg, " is not installed. Exiting...\n")); quit(status = 1)} else {cat(paste0(pkg, " is installed.\n"))}}'; then
    exit 1
fi

echo -e "\nChecking python dependencies ...\n"
# Check Python modules
if ! python -c 'import sys; modules = ["sys", "pandas", "os"]; [print(f"{mod} is installed.") if __import__(mod) else (print(f"{mod} is not installed. Exiting...\n"), sys.exit(1)) for mod in modules]'; then
    exit 1
fi

echo -e "\nAll dependencies are installed.\n"



echo -e "\nStarting prediction ...\n"
echo -e "\nStarting prediction ...\n" >> "$output"/program.log

# Launch the predictions from every fasta and fna files stored in the input folder
for file in "${input_folder[@]}"; do 
    # Determine the additional arguments based on the presence of -H and -D
    database_arg="$database_file"
    hmm_arg=0
    nameFile=$(basename $file)

    if [ ! -z "$database_secondary" ]; then
        database_arg="$database_secondary"  # Use secondary database if specified
    fi

    if [ ! -z "$perfil_HMM" ]; then
        hmm_arg="$perfil_HMM"  # Use HMM profile if specified
    fi

    # Execute prediction script with the appropriate arguments
    $SCRIPT_DIR/prediction.sh "$file" "$database_arg" "$hmm_arg" "$output" >> "$output"/program.log 2>&1

    # Check if the prediction was successful
    if [ $? -ne 0 ]; then
      echo "Prediction failed for file $file." >> "$output"/program.log
      echo -e "{$RED}Prediction failed{$NC} for $nameFile. Check logs for details."
    else    
      echo -e "${GREEN}Prediction completed${NC} for $nameFile.\n"
    fi
done

## Formatting the results: output of the prediction.sh is: database/total_proteins.fasta
### the original header is: >GCA_002245475_Pararge_aegeria_120629693_basic_juvenile_hormone-suppressible_protein_1-like|NJDC01000609.1:143746-344643|-|1212|0|5|195946|200866|200866[200866]:#200768[200768]:99[99]|200018[200018]:199884[199884]:135[135]|198727[198727]:197900[197900]:828[828]|197622[197622]:197449[197449]:174[174]|196971[196959]:195946[195946]:1026[1014]

echo -e "\nStarting formatting the results ...\n" 
echo -e "\nStarting formatting the results ...\n" >> "$output"/program.log

[ -s "$output/total_proteins.fasta" ] || { echo -e "\nPrediction file is empty. Exiting ..."; exit 1; }


# The header will be after awk:  >GCA_002245475_NJDC01000609_143746-344643_-_0_195946_200866 (after coordinates we have: strand, e-value, start and stop)
awk -F '|' '{if ($1 ~ /^>GCA/) print $1"_"$2"_"$3"_"$5"_"$7"_"$8; else print $0}' "$output"/total_proteins.fasta | sed 's/\..:/_/g' | awk -F '_' '{if ($1 ~ /^>GCA/) print $1"_"$2"_"$(NF-5)"_"$(NF-4)"_"$(NF-3)"_"$(NF-2)"_"$(NF-1)"_"$NF; else print $0}' > "$output"/total_proteins_header.fasta

# Remove duplicate by name
seqkit rmdup -n "$output"/total_proteins_header.fasta > "$output"/proteins_rmdup.fasta
#remove temporary file 

[ -f "$output"/total_proteins_header.fasta ] && rm "$output"/total_proteins_header.fasta

# Remove protein longer and shorter than x (it will be up to the client choose this parameter) -m -M
if [ ! -z "$min_length" ] || [ ! -z "$max_length" ]; then
    echo -e "\n Filtering proteins with ${min_length:+Min_length $min_length} ${max_length:+Max_length $max_length} ...\n" 
    echo -e "\n Filtering proteins with ${min_length:+Min_length $min_length} ${max_length:+Max_length $max_length} ...\n" >> "$output"/program.log
    seqkit seq ${min_length:+-m $min_length} ${max_length:+-M $max_length} "$output"/proteins_rmdup.fasta > "$output"/total_protein_header_filtered.fasta
    [ -s "$output"/total_protein_header_filtered.fasta ] || { echo -e "\nNo proteins after filtering by size. Exiting ..."; exit 1; }
    [ -f  "$output"/proteins_rmdup.fasta ] && rm "$output"/proteins_rmdup.fasta  # Remove intermediate file
else
    mv "$output"/proteins_rmdup.fasta "$output"/total_protein_header_filtered.fasta
fi


# Filter by domain (up to user). IMPORTANT ADD INTERPROSCAN IN .bashrc
if [[ -n "$domain" ]]; 
then
  echo -e "\nRunning interproscan on $domain domain(s) ...\n" 
  echo -e "\nRunning interproscan on $domain domain(s) ...\n" >> "$output"/program.log
  interproscan.sh -i "$output"/total_protein_header_filtered.fasta -b "$output"/interpro_results_database -f tsv  >> "$output"/program.log 2>&1
  [ -d $PWD/temp ] && rm -r $PWD/temp/
  #Python script to keep only the protein with the specific domain in th -D argument
  python3 $SCRIPT_DIR/find_protein_interpro.py "$output"/interpro_results_database.tsv "$output"/proteins_to_filter.txt "$domain" 
  [ -s "$output"/proteins_to_filter.txt ] ||  { echo -e "\nNo proteins after filtering by domain. Exiting ..."; exit 1; }
  # Take only the filtered protein by name
  seqkit grep -f "$output"/proteins_to_filter.txt "$output"/total_protein_header_filtered.fasta -o "$output"/filtered_protein.fasta
  [ -f  "$output"/interpro_results_database.tsv ] && rm "$output"/interpro_results_database.tsv	
  [ -f  "$output"/interpro_results_database_1.tsv ] && rm "$output"/interpro_results_database_1.tsv	
  [ -f  "$output"/proteins_to_filter.txt ] && rm "$output"/proteins_to_filter.txt
else
  mv "$output"/total_protein_header_filtered.fasta "$output"/filtered_protein.fasta
fi


# amminoacid sequences in one line in order to save the ouput such a .csv file to formatting it in R. The output will be ordered and analyzed in R (delete overlap region)
awk -F '_' '{if ($1 ~ /^>GCA/) print $1"_"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8; else print $0}'  "$output"/filtered_protein.fasta | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' > "$output"/filtered_protein.tsv

# Launch script in R in order to take off all the overlapping region
# The output of the script is a txt file with the name of the best sequences for region and a file excel with the results
Rscript $SCRIPT_DIR/order_data.R "$output"

# Take off e-value for seqkit grep
# Take off e-value for seqkit grep from the fasta file
awk -F '_' '{if ($1 ~ /^>GCA/) print $1"_"$2"_"$3"_"$4"_"$5"_"$7"_"$8; else print $0}' "$output"/filtered_protein.fasta > "$output"/protein_to_filter.fasta
# Take off e-value. This step is usefull for seqkit (grep) from the txt generated with R 
awk -F '_' '{if ($1 ~ /GCA/) print $1"_"$2"_"$3"_"$4"_"$5"_"$7"_"$8; else print $0}' "$output"/filtered_protein.txt >  "$output"/filtered_protein1.txt

seqkit grep -n -f  "$output"/filtered_protein1.txt  "$output"/protein_to_filter.fasta > "$output"/filtered_protein_final.fasta

####   Phylogenic reconstruction
if [[ $iqtree = true ]]; 
then 
  mkdir -p "$output"/alignment
  sed -e 's/-/_/g;s/_+_/_/g;s/___/_/g' "$output"/filtered_protein_final.fasta > "$output"/alignment/final_protein_species_alignment.fasta
  num_prot=$(grep -c ">" "$output"/alignment/final_protein_species_alignment.fasta)
  if [[ $num_prot > 1 ]];
  then
    echo -e "\nRunning kalign on $num_prot proteins ...\n"
    echo -e "\nRunning kalign on $num_prot proteins ...\n" >> "$output"/program.log
    kalign -i "$output"/alignment/final_protein_species_alignment.fasta --format fasta -o "$output"/alignment/kalign.fasta
    echo -e "\nRunning trimmal ...\n"
    echo -e "\nRunning trimmal ...\n" >> "$output"/program.log
    trimal -in "$output"/alignment/kalign.fasta -fasta -out "$output"/alignment/clean_kalign.fasta -automated1
    echo -e "\nRunning iqtree ...\n"
    echo -e "\nRunning iqtree ...\n" >> "$output"/program.log
    $SCRIPT_DIR/iqtree.sh -s "$output"/alignment/clean_kalign.fasta $iq_params >> "$output"/program.log 2>&1
    [ -f "$output"/alignment/final_protein_species_alignment.fasta ] && rm "$output"/alignment/final_protein_species_alignment.fasta
    [ -f "$output"/alignment/kalign.fasta ] && rm "$output"/alignment/kalign.fasta
  else
    echo -e "\nCannot run kalign on 1 protein! Skipping ...\n"
  fi
fi

# Retrieve species name if the file was downloaded from NCBI
if [[ "$species_name" = true ]]; 
then 
  $SCRIPT_DIR/species_name.sh "$output"/filtered_protein_final.fasta "$output"/final_protein_species.fasta
  ### save the result as csv
  awk -F '_' '{if ($1 ~ /^>GCA/) print $1"_"$2","$3","$4","$5","$6","$7","$8"_"$9}'  "$output"/final_protein_species.fasta  > "$output"/final_protein_species.csv
  Rscript $SCRIPT_DIR/save_result_species.R "$output"
  [ -f  "$output"/final_protein.xlsx ] && rm "$output"/final_protein.xlsx
  [ -f  "$output"/filtered_protein.tsv ] && rm "$output"/filtered_protein.tsv
  [ -f  "$output"/final_protein_species.csv ] && rm "$output"/final_protein_species.csv
fi
# Remove files
[ -f "$output"/filtered_protein.tsv ] && rm "$output"/filtered_protein.tsv
[ -f "$output"/filtered_protein1.txt ] && rm "$output"/filtered_protein1.txt
[ -f "$output"/filtered_protein.fasta ] && rm "$output"/filtered_protein.fasta
[ -f "$output"/protein_to_filter.fasta ] && rm "$output"/protein_to_filter.fasta

echo -e "\n\n${GREEN}SUCCESS${NC} - Pipeline completed."