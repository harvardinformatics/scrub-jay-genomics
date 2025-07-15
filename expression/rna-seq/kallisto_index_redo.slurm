#!/bin/bash
#SBATCH -n 2                # Number of cores
#SBATCH -N 1
#SBATCH -t 3-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared,holy-smokes   # Partition to submit to
#SBATCH --mem=32G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o index.out # File to which STDOUT will be written
#SBATCH -e index.err # File to which STDERR will be written

source ~/mambaforge/bin/activate kallisto

cd /n/holylfs05/LABS/informatics/Lab/scrubjay/GENOMIC_RESOURCES/gene_annotation/results/fixed_merged_TOGA

kallisto index -i transcripts_merged_TOGA_tgutRef_noDup.kalIdx transcripts_merged_TOGA_tgutRef_noDup.fa
