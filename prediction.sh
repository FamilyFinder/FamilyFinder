#!/bin/bash

#SBATCH -n 20
#SBATCH -p long
#SBATCH -t 5-20:00:00
#SBATCH -J anotación
#SBATCH -e error.err
#SBATCH -o output_miniprot

################## inputs ####################
input=$1
proteins_database=$2
database_protein_1=$3
output=$4
hmm_file=$5
home=$HOME
genome_dir=$(dirname $input)

# Take the code of the genome. Example GCA_0001234
species_name=$(basename "$input" | sed -e 's/.fasta//g;s/.fna//g' | awk -F '.' '{print $1}')

echo -e '1. Analyzing' $species_name 'with miniprot ...' >> "$output"/program.log
mkdir -p "$output"/prediction_"$species_name"/miniprot
mkdir -p "$output"/prediction_"$species_name"/metaeuk/orthodb
mkdir -p "$output"/prediction_"$species_name"/metaeuk/orthodb1
mkdir -p "$output"/prediction_"$species_name"/metaeuk/orthodb2
mkdir -p "$output"/prediction_"$species_name"/metaeuk/hmmscan

# If there are two databases as inputs make another directory where to save the results
if [ $database_protein_1 != 0 ]; then 
    mkdir -p "$output"/prediction_"$species_name"/metaeuk/databases_1; 
fi
cd "$output"/prediction_"$species_name"/miniprot



################ 1. miniprot ################
#############################################

# miniprot -t8 -d output.mpi genoma_da_index.fna
echo "Running miniprot on $input ..."  >> "$output"/program.log
[[ -f "$species_name".mpi ]] || miniprot -t 16 -d "$species_name".mpi "$input"
# $2= reference proteins database
miniprot -G 600k --gff "$species_name".mpi "$proteins_database" > "$species_name".gff


################ 1.1 Formatting output ################
#######################################################

echo -e '1.1 Formatting miniprot output ...\n' >> "$output"/program.log
# Sorting by coordinates find with miniprot
grep 'mRNA'  "$species_name".gff | sort -k1,1 -k4,4n -k5,5n > "$species_name"_filtered.gff
# Make only one big region taking off the overlapping regions
bedtools merge -i "$species_name"_filtered.gff > "$species_name"_filtered.bed
# IMPORTANT: Change and check and put $3  += 7000 (if $3 > max of the contig)
awk 'BEGIN{OFS="\t"} {$2 -= 7000} {print $0}' "$species_name"_filtered.bed | awk -F '\t' '{if($2<0) print $1 "\t" "1" "\t" $3; else print $0}' > "$species_name"_final_filtered.bed
# Taking only the region of interest
bedtools getfasta -fi "$input" -bed "$species_name"_final_filtered.bed -fo "$output"/prediction_"$species_name"/"$species_name"_filtered_NoDup.fasta
cd ..

################ 2. Metaeuk prediction ################
#######################################################

echo -e '2. Starting Metaeuk prediction ...\n' >> "$output"/program.log

# if there are two databases as input make this additional prediction
if [ $database_protein_1 != 0 ]; then 
    metaeuk easy-predict "$species_name"_filtered_NoDup.fasta "$database_protein_1" metaeuk/databases_1/prediction_"$species_name"_database1 metaeuk/databases_1 \
    --threads 20 --remove-tmp-files 1 --max-intron 600000 --compressed 1; 
fi

# prediction with orthodb
metaeuk easy-predict "$species_name"_filtered_NoDup.fasta "$proteins_database" metaeuk/orthodb/prediction_"$species_name"_ortho metaeuk/orthodb \
    --threads 20 --remove-tmp-files 1 --max-intron 600000  -s 7.5 

# make another prediction with orthodb because changing the sensibility we can find DIFFERENTS proteins (little diferrences): ORTHODB1
metaeuk easy-predict "$species_name"_filtered_NoDup.fasta "$proteins_database" metaeuk/orthodb1/prediction_"$species_name"_ortho1 metaeuk/orthodb1 \
--threads 20 --remove-tmp-files 1 --max-intron 600000 

# make another prediction with orthodb because changing the sensibility we can find DIFFERENTS protein (little diferrences): ORTHODB2
metaeuk easy-predict "$species_name"_filtered_NoDup.fasta "$proteins_database" metaeuk/orthodb2/prediction_"$species_name"_ortho2 metaeuk/orthodb2 \
--threads 20 --remove-tmp-files 1 --max-intron 600000 -s 1.0

# copy the usefull files
mv metaeuk/*/prediction_"$species_name"_*.gff "$output"/prediction_"$species_name"/



################ 3. Formatting Metaeuk outputs ################
###############################################################
# Filter by taking uniq entry. The columns select are the start and stop column.

echo -e 'Formatting Metaeuk outputs ...\n'  >> "$output"/program.log

databases=("orthodb" "orthodb1" "orthodb2")

# Loop through each db
for db in "${databases[@]}"; do
    # Perform deduplication
    seqkit rmdup -s "metaeuk/$db/prediction_${species_name}_ortho${db#orthodb}.fas" > "metaeuk/$db/prediction_${species_name}_filtered_ortho${db#orthodb}.fasta"
        
    # Extract and sort unique protein identifiers
    grep '>' "metaeuk/$db/prediction_${species_name}_filtered_ortho${db#orthodb}.fasta" | sort -t '|' -k2,2 -k7,7 -u | sort -t '|' -k1,1 -u | sed 's/>//g' > "proteins_to_filter${db#orthodb}.txt"
        
    # Filter the fasta file for unique proteins
    seqkit grep -f "proteins_to_filter${db#orthodb}.txt" "metaeuk/$db/prediction_${species_name}_filtered_ortho${db#orthodb}.fasta" -n >> proteins_orthodb_metaeuk.fasta
done

# if there is a second databases make the filter as the other. 
if [ -d "metaeuk/databases_1" ]; then
    seqkit rmdup -s "metaeuk/databases_1/prediction_${species_name}_database1.fas" > "metaeuk/databases_1/prediction_${species_name}_database1_filtered.fasta"
    grep '>' "metaeuk/databases_1/prediction_${species_name}_database1_filtered.fasta" | sort -t '|' -k2,2 -k7,7 -u | sort -t '|' -k1,1 -u | sed 's/>//g' > "proteins_to_filter_database1.txt"
    seqkit grep -f "proteins_to_filter_database1.txt" "metaeuk/databases_1/prediction_${species_name}_database1_filtered.fasta" -n >> proteins_orthodb_metaeuk.fasta
else
    echo -e "No optional database found.\n"  >> "$output"/program.log
fi

# if the option hmm is active. filter by hmm 
if [[ ! -z "$hmm_file" ]]; 
then 
    hmmscan --tblout metaeuk/hmmscan/hmmscan "$hmm_file" proteins_orthodb_metaeuk.fasta 
    # Filter
    awk '{print $3}' metaeuk/hmmscan/hmmscan  | head -n -10 | tail -n +4 > metaeuk/hmmscan/hmm_filtered.txt
    # Taking unique entry 
    sort -t '|' -k1,1 -u metaeuk/hmmscan/hmm_filtered.txt | sort -t '|' -k2,2 -k7,7 -u  > metaeuk/hmmscan/proteins_to_filter_hmm.txt 
    seqkit grep -f metaeuk/hmmscan/proteins_to_filter_hmm.txt proteins_orthodb_metaeuk.fasta -n > proteins_orthodb_metaeuk_filtered.fasta
    mv metaeuk/hmmscan/hmmscan hmmscan
else
    cp proteins_orthodb_metaeuk.fasta proteins_orthodb_metaeuk_filtered.fasta
fi

############### Formatting the output with all the protein from all the predictions ###############
sed -i "s/>/&$species_name\_/" proteins_orthodb_metaeuk_filtered.fasta
# remove duplicated (based on sequences)
seqkit rmdup -s proteins_orthodb_metaeuk_filtered.fasta > filtered_proteins.fasta
# Remove protein with 90% of identity
cd-hit -i filtered_proteins.fasta  -o final_proteins_"$species_name".fasta -n 5 -c 0.9 -aL 0.9 

# Copy the protein in a uniq file to make a phylogenetic reconstruction

cat final_proteins_"$species_name".fasta >> "$output"/total_proteins.fasta
## Save nucleotidic sequences
mv "$species_name"_filtered_NoDup.fasta "$species_name"_nucleotidic_sequences.fasta 
# Remove all the useless files
[ -f proteins_to_filter_1.txt ] && rm proteins_to_filter_1.txt
[ -f proteins_to_filter_2.txt ] && rm proteins_to_filter_2.txt
if [ $database_protein_1 != 0 ]; then 
    rm proteins_to_filter_database1.txt; 
fi 
rm -r metaeuk/
rm -r miniprot/
[ -f filtered_proteins.fasta ] && rm filtered_proteins.fasta
[ -f proteins_orthodb_metaeuk.fasta ] && rm proteins_orthodb_metaeuk.fasta
[ -f proteins_orthodb_metaeuk_filtered.fasta ] && rm proteins_orthodb_metaeuk_filtered.fasta
[ -f final_proteins_"$species_name".fasta.clstr ] && rm final_proteins_"$species_name".fasta.clstr


