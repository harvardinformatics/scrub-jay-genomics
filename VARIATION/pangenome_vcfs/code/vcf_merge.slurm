#!/bin/bash
#SBATCH -J sort # A single job name for the array
#SBATCH --partition=holy-smokes,shared # Partition
#SBATCH -n 4 # Number of Nodes required
#SBATCH --mem=32G # Memory request
#SBATCH -t 3-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o sort_%J.out # Standard output
#SBATCH -e sort_%J.err # Standard error


source ~/mambaforge/bin/activate bcftools

bcftools merge --force-samples -O z -o ../results/pgg_MERGED.vcf.gz ../results/*bub_wave_sort.vcf.gz
