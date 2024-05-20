#!/bin/bash

#SBATCH -n 30
#SBATCH --mem=20000
#SBATCH -p long
#SBATCH -t 10-20:00:00
#SBATCH -J anotación
#SBATCH -e error.err
#SBATCH -o output


iqtree "$@" -m TEST 
