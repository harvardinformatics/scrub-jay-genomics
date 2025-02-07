#!/bin/bash
#SBATCH -J sort # A single job name for the array
#SBATCH -a 1-44
#SBATCH --partition=holy-cow,shared # Partition
#SBATCH -n 4 # Number of Nodes required
#SBATCH --mem=16G # Memory request
#SBATCH -t 3-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o sort_%J.out # Standard output
#SBATCH -e sort_%J.err # Standard error


VCF=$(sed "${SLURM_ARRAY_TASK_ID}q;d" vcf_fofn.txt)
DIR="$(dirname "${VCF}")"
SAMPLE="$(basename -s .vcf "${VCF}")"

cd $DIR

source ~/mambaforge/bin/activate bcftools

#bcftools sort -m 16G -o ${SAMPLE}_sort.vcf.gz -O z $VCF
bcftools index ${SAMPLE}_sort.vcf.gz
