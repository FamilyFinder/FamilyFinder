### FamilyFinder: a new pipeline for the rapid retrieval of a protein or a protein family in an unannotated genomes
If you use this pipeline please cite:

Familyfinder takes an assembled genome as input to accurately predict a set or protein of interest. This pipeline uses several other well-known tools, such as miniport and metaeuk, to extract only the nucleotide sequences of interest from a given genome. In addition, when the pipeline was tested on insect genomes, it took about one minute to extract the protein and nucleotide sequences from the genome. 

Usage: Family_finder.sh -i <input> -o <output> -d <database> [options]

The mandatory Arguments are:
-i <input>         Path to the directory containing the genomes.
-d <database>      Database for protein predictions (e.g., orthoDB). This database
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



In order to use the pipeline you need to have these dependencies installed:
-BUSCO=5.7.1 
-HMMer (if you want to filter out the protein thanks a personal profile HMM. This it's up to user to make it) 
-miniprot
-seqkit
-bedtools
-iqtree (if the option -t is setted)
-entrez Direct (if the option -s in setted)
-R (within packages dplyr,tidyr, writexl, stringr and data.table)
-Python v 3.x
-Interproscan (if the option -I is setted)


<!--
**FamilyFinder/FamilyFinder** is a ✨ _special_ ✨ repository because its `README.md` (this file) appears on your GitHub profile.

Here are some ideas to get you started:

- 🔭 I’m currently working on ...
- 🌱 I’m currently learning ...
- 👯 I’m looking to collaborate on ...
- 🤔 I’m looking for help with ...
- 💬 Ask me about ...
- 📫 How to reach me: ...
- 😄 Pronouns: ...
- ⚡ Fun fact: ...
-->
