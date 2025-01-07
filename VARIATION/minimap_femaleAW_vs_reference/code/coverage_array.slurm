#!/bin/bash
#SBATCH -J coverage # A single job name for the array
#SBATCH -a 1-7
#SBATCH --partition=shared # Partition
#SBATCH -n 2 # Number of Nodes required
#SBATCH -N 1
#SBATCH --mem=48G # Memory request 
#SBATCH -t 3-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o coverage_%J_%a.out # Standard output
#SBATCH -e coverage_%J_%a.err # Standard error

source ~/mambaforge/bin/activate bedtools


#Assign assembly and read files from spec sheet
BAM=$(sed "${SLURM_ARRAY_TASK_ID}q;d" bam_fofn.txt)
ID=$(basename $BAM .mapped.sorted.bam)

bedtools coverage -sorted -hist -b $BAM -a aphWoo1_1kb_windows_chZ.bed > ${ID}_chZ_1kb_cov.txt
