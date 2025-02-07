#!/bin/bash
#SBATCH -n 2                # Number of cores
#SBATCH -t 3-0:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -a 1-45
#SBATCH -p holy-smokes   # Partition to submit to
#SBATCH --mem=32G           # Memory pool for all cores (see also --mem-per-cpu)

source ~/mambaforge/bin/activate bcftools

sample=$(sed "${SLURM_ARRAY_TASK_ID}q;d" pggb_samples.txt)

bcftools view -Oz -s $sample -o SPLIT_pggb/${sample}_pggb.vcf.gz pggb_merged.vcf

#for file in svim_merged.nameFix.vcf.gz; do
#  for sample in `bcftools query -l $file`; do
#    bcftools view  -Oz -s $sample -o SPLIT_svim/${sample}_svim.vcf.gz $file
#  done
#done

#for file in sj-88b.call.nameFix.reHeader.vcf.gz; do
#  for sample in `bcftools query -l $file`; do
#    bcftools view  -Oz -s $sample -o SPLIT_minigraph/${sample}_minigraph.vcf.gz $file
#  done
#done

