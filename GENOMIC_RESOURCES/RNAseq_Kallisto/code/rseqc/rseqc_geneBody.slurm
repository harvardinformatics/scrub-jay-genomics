#!/bin/bash
#SBATCH -J rseqc # A single job name for the array
#SBATCH --partition=holy-cow # Partition
#SBATCH -n 2 # Number of Nodes required
#SBATCH -N 1
#SBATCH --mem=48G # Memory request 
#SBATCH -t 14-00:00 # Maximum execution time (D-HH:MM)
#SBATCH -o rseqc.out # Standard output
#SBATCH -e rseqc.err # Standard error

source ~/mambaforge/bin/activate rseqc

geneBody_coverage.py -r toga_temp.bed -i bam_list.txt -o geneBody_coverage
