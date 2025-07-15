#!/bin/bash
#SBATCH -J vg_decon # A single job name for the array
#SBATCH -a 1-46
#SBATCH --partition=shared # Partition
#SBATCH -n 4 # Number of Nodes required
#SBATCH --mem=32G # Memory request 
#SBATCH -t 7-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o vgde_%A_%a.out # Standard output
#SBATCH -e vgde_%A_%a.err # Standard error


COM=$(sed "${SLURM_ARRAY_TASK_ID}q;d" community_gfa_fofn.txt)
OUTDIR=$(echo $COM | awk 'BEGIN{FS="/"; OFS="/"} {print $1,$2,$3}')
SAMPLE=$(basename $COM .gfa)

sed 's/AW_365336#REF/aphWoo1#REF/g' $COM > ${OUTDIR}/${SAMPLE}_nameFix.gfa

./vg deconstruct -P aphWoo1 -H '#' -e -a -t 4 ${OUTDIR}/${SAMPLE}_nameFix.gfa  > ${OUTDIR}/${SAMPLE}_nameFix.vcf
