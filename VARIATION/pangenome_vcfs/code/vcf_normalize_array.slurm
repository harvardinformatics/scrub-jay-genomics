#!/bin/bash
#SBATCH -J vg_decon # A single job name for the array
#SBATCH -a 1-46
#SBATCH --partition=shared # Partition
#SBATCH -n 8 # Number of Nodes required
#SBATCH --mem=64G # Memory request 
#SBATCH -t 7-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o vcfnorm_%A_%a.out # Standard output
#SBATCH -e vcfnorm_%A_%a.err # Standard error


COM=$(sed "${SLURM_ARRAY_TASK_ID}q;d" community_vcf_fofn.txt)
OUTDIR=$(echo $COM | awk 'BEGIN{FS="/"; OFS="/"} {print $1,$2,$3}')
SAMPLE=$(basename $COM .vcf)

gzip $COM

singularity exec pggb_native/pggb_cannon.sif vcfbub -l 0 -a 100000 --input ${COM}.gz > ${OUTDIR}/${SAMPLE}_nameFix_bub.vcf
singularity exec pggb_native/pggb_cannon.sif vcfwave -I 1000 -t 8 ${OUTDIR}/${SAMPLE}_nameFix_bub.vcf > ${OUTDIR}/${SAMPLE}_nameFix_bub_wave.vcf
