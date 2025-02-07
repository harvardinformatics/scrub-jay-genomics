#!/bin/bash
#SBATCH -n 1                # Number of cores
#SBATCH -t 3-0:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared   # Partition to submit to
#SBATCH --mem=24G           # Memory pool for all cores (see also --mem-per-cpu)


source ~/mambaforge/bin/activate bcftools

bcftools view -H --types indels -i 'LEN>=50' ../results/pggb_merged.vcf | awk 'BEGIN{FS="\t"; OFS="\n"} length($4)<length($5){print ">"$1":"$2"_"NR,$5} length($4)>length($5){print ">"$1":"$2"_"NR,$4}' > ../results/pggb_merged_INDEL_alleles.fasta
