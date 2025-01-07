#!/bin/bash
#SBATCH -J mapping # A single job name for the array
#SBATCH --partition=test # Partition
#SBATCH -n 2 # Number of Nodes required
#SBATCH -N 1
#SBATCH --mem=24G # Memory request 
#SBATCH -t 0-08:00 # Maximum execution time (D-HH:MM)
#SBATCH -o mapping.out # Standard output
#SBATCH -e mapping.err # Standard error

source ~/mambaforge/bin/activate samtools

for bam in ../results/stringtie_RNAseq/WORKFLOW/*_hisat/*.bam
do
	echo $bam 
	samtools flagstat $bam
done
