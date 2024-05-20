#!/bin/bash


#SBATCH -n 20
#SBATCH -p short
#SBATCH -t 00-12:00:00
#SBATCH -J anotación
#SBATCH -e /mnt/lustre/home/fblasio/error.err
#SBATCH -o /mnt/lustre/home/fblasio/output_miniprot

# Entrez tools. This command allow to take the species name for a specific entry (example GCA_0001..)
while read -u 9 line;
do 
	if [[ "$line" == ">"GCA* ]] ;
	then
		assembly=$(echo "$line" | awk -F '_' '{print $1"_"$2}' | sed 's/>//g')
		name=$(esearch -db assembly -query "$assembly" | efetch -format docsum | xtract -pattern DocumentSummary -element SpeciesName | uniq)
		echo "$line""_""$name"
	else 
		echo "$line"
	fi
done 9< $1 > $2

sed -i "s/ /_/g" $2
