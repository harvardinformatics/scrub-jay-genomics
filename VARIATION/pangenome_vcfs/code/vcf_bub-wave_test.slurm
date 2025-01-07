#!/bin/bash
#SBATCH -J vg_decon # A single job name for the array
#SBATCH --partition=holy-cow # Partition
#SBATCH -n 8 # Number of Nodes required
#SBATCH --mem=64G # Memory request 
#SBATCH -t 7-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o vcfnorm_test.out # Standard output
#SBATCH -e vcfnorm_test.err # Standard error


VCF="../results/allbird_community.1.vcf.gz"

singularity exec pggb_11-15-23.sif vcfbub -l 0 -a 100000 --input $VCF > test_bub.vcf
singularity exec pggb_11-15-23.sif vcfwave -I 1000 -t 8 test_bub.vcf > test_bub-wave.vcf
