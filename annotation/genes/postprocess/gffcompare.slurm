#!/bin/bash
#SBATCH -J gffcompare # A single job name for the array
#SBATCH --partition=shared # Partition
#SBATCH -n 2 # Number of Nodes required
#SBATCH -N 1
#SBATCH --mem=42G # Memory request
#SBATCH -t 3-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o gffcompare.out # Standard output
#SBATCH -e gffcompare.err # Standard error

TGUT="scrubjay_tGut_regID.gtf"
GGAL="scrubjay_gGal_regID.gtf"

/n/holylfs05/LABS/informatics/Users/dkhost/gffcompare-0.12.6.Linux_x86_64/gffcompare -o merged_TOGA_noDup -D -r $TGUT $GGAL

