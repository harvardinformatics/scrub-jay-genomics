#!/bin/bash
#SBATCH -J quast_array # A single job name for the array
#SBATCH -a 1-11
#SBATCH --partition=holy-smokes # Partition
#SBATCH -n 4 # Number of Nodes required
#SBATCH --mem=24G # Memory request 
#SBATCH -t 1-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o quast_%A_%a.out # Standard output
#SBATCH -e quast_%A_%a.err # Standard error

source activate assembly_eval

p=$(sed "${SLURM_ARRAY_TASK_ID}q;d" assemblies_fofn.txt)

sample=$(basename $p .fa.gz)

quast.py -t 4 --eukaryote -o ${sample}_quast $p
