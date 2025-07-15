#!/bin/bash

#SBATCH -N 1
#SBATCH -t 14-00:00:00
#SBATCH --mem 400G
#SBATCH -n 8
#SBATCH --partition holy-cow
#SBATCH -J snakemake

source ~/mambaforge/bin/activate snakemake

snakemake -r --cores 8 --use-conda --rerun-incomplete
