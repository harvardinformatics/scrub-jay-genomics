#!/bin/bash
#SBATCH -n 2                # Number of cores
#SBATCH -t 3-0:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -a 1-45
#SBATCH -p holy-cow   # Partition to submit to
#SBATCH --mem=32G           # Memory pool for all cores (see also --mem-per-cpu)


source ~/mambaforge/bin/activate bcftools

SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" samples.txt)

bcftools view -o MERGED/${SAMPLE}_TEMP_pggb.vcf SPLIT_pggb/${SAMPLE}_pggb.vcf.gz

bcftools view -o MERGED/${SAMPLE}_TEMP_svim.vcf SPLIT_svim/${SAMPLE}_svim.vcf.gz

bcftools view -o MERGED/${SAMPLE}_TEMP_mini1.vcf SPLIT_minigraph/${SAMPLE}#1_minigraph.vcf.gz

bcftools view -o MERGED/${SAMPLE}_TEMP_mini2.vcf SPLIT_minigraph/${SAMPLE}#2_minigraph.vcf.gz

ls MERGED/${SAMPLE}_TEMP*.vcf > ${SAMPLE}_vcf_fofn.txt

/n/holylfs05/LABS/informatics/Users/dkhost/SURVIVOR-master/Debug/SURVIVOR merge ${SAMPLE}_vcf_fofn.txt 1000 1 1 1 0 50 MERGED/${SAMPLE}_merged.vcf

rm MERGED/${SAMPLE}_TEMP*.vcf
