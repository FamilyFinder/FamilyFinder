# FamilyFinder: a new pipeline for the rapid retrieval of a protein or a protein family in an unannotated genomes
If you use this pipeline please cite: XXXXXXX

Familyfinder takes an assembled genome as input to accurately predict a set or protein of interest. This pipeline uses several other well-known tools, such as miniprot and metaeuk, to extract only the nucleotide sequences of interest from a given genome. In addition, when the pipeline was tested on insect genomes, it took about one minute to extract the protein and nucleotide sequences from the selected genome. 

The figure below shows the workflow of the pipeline.
<br />  <img src=https://github.com/FamilyFinder/FamilyFinder/assets/170311637/71728427-d9ea-42da-9f9b-1f85c7a5c37b width="600" height="700">






<br /> Usage: Family_finder.sh ``-i <input> -o <output> -d <database> [options]``

**The mandatory Arguments are**:
- -i ``<input>``; Path to the directory containing the genomes.
- -d ``<database>``; Database for protein predictions (e.g., orthoDB). This database will be used with miniprot (https://github.com/lh3/miniprot).
  
**The optional Arguments are**:
- -o ``<output>``; Output directory (default: ./Family_finder_out).
- -D ``<database>``; Secondary database for predictions, not used with miniprot. Used alongside -d database by metaeuk.
- -m ``<min_length>``; Minimum protein length to consider after predictions.
- -M ``<max_length>``; Maximum protein length to consider after predictions.
- -H ``<file.hmm>``; HMM file to filter prediction results (recommended when using -D).
- -t  Specify to generate a phylogenetic tree after predictions with IQtree. 
- -q ``<iqtree_params>``; Specify iqtree parameters (default: -bb 1000 -alrt 1000 -abayes -mem 20G -safe -T AUTO; http://www.iqtree.org/doc/).
- -s  Specify to retrieve species names. **Note**: This may significantly increase processing time.
- -I ``<domain>``; Specify domains using InterPro classification, separated by commas (e.g., IPR001,IPR002,IPR004).
- -e ``<environment>``; Specify a Conda environment for third-party components.
- -h  Display this help message and exit.

**Examples**: 
<br /> BASIC: ``Family_finder.sh -i /path/to/genomes -d path/to/orthoDB -o /path/to/output``
<br /> WITH PHYLOGENETIC RECONSTRUCTION: ``Family_finder.sh -i /path/to/genomes -d path/to/orthoDB -o /path/to/output -t``
<br /> WITH FILTER BY PROTEIN DOMAIN: ``Family_finder.sh -i /path/to/genomes -d path/to/orthoDB -o /path/to/output -I IPR001,IPR002,IPR004``

**Note:** There is a special script to filter the protein by domain and size when you download the protein from othoDB. The script name is Formatting_database.sh. Due to the large header in orthoDB, it is recommended to use this script if you download the protein from this database.
<br /> The usage of the script is ``Formatting_database.sh -i <input> -o <output> [options]``
<br /> **Mandatory arguments**:
- -i ``<input>``; Input file in FASTA format. The first line must start with '>'.
- -o ``<output>``; Output file.
  
**Optional Arguments:**
- -m ``<min_length>``; Minimum protein length to consider in the database.
- -M ``<max_length>``; Maximum protein length to consider in the database.
- -I ``<domain>``; Filter by domain. Specify domains using InterPro classification, separated by commas (e.g., IPR001,IPR002,IPR004).
- -h ; Display this help message and exit.

**Examples**:
<br /> ``Formatting_databases.sh -i input.fasta -o output.fasta``
<br /> ``Formatting_databases.sh -i input.fasta -o output.fasta -m 200 -M 1000 -I IPR001,IPR002,IPR004``

You need to have these dependencies installed in order to use the pipeline (The versions refers to the version that were used for the test of the pipeline.):
- BUSCO=5.7.1
- Metaeuk v. 6.a5d39d9
- HMMer 3.x (if you want to filter out the protein thanks a personal profile HMM. This it's up to user to make it) 
- Miniprot v. 0.7
- Seqkit
- Bedtools
- Iqtree (if the option -t is selected)
- Entrez Direct (if the option -s in selected)
- R (within packages dplyr,tidyr, writexl, stringr and data.table)
- Python v 3.x
- Interproscan v.  5.65-97.0 (if the option -I is selected)

To test the pipeline please download the test_data folder from https://drive.google.com/drive/u/1/folders/1f7NsnSioCTsnZ0Uh978qSA-_C506eQka

