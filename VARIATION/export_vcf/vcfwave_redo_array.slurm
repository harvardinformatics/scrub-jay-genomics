#!/bin/bash
#SBATCH -J vcfwave # A single job name for the array
#SBATCH -a 1-46
#SBATCH --partition=shared # Partition
#SBATCH -n 8 # Number of Nodes required
#SBATCH --mem=64G # Memory request 
#SBATCH -t 3-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o vcfwavw_%A_%a.out # Standard output
#SBATCH -e vcfwave_%A_%a.err # Standard error

source ~/mambaforge/bin/activate vcflib-manual

COM=$(sed "${SLURM_ARRAY_TASK_ID}q;d" community_vcf_fofn.txt)
OUTDIR=$(echo $COM | awk 'BEGIN{FS="/"; OFS="/"} {print $1,$2,$3}')
SAMPLE=$(basename $COM .vcf)

#gzip $COM
mkdir $OUTDIR/vcfwave_redo

#singularity exec pggb_native/pggb_cannon.sif vcfbub -l 0 -a 100000 --input ${COM}.gz > ${OUTDIR}/${SAMPLE}_nameFix_bub.vcf
/n/holylfs05/LABS/informatics/Users/dkhost/vcflib/build/vcfwave -I 1000 -t 8 ${OUTDIR}/${SAMPLE}_nameFix_bub.vcf > ${OUTDIR}/vcfwave_redo/${SAMPLE}_nameFix_bub_wave.vcf
