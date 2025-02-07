#!/bin/bash
#SBATCH -J miniamp # A single job name for the array
#SBATCH -a 1-7
#SBATCH --partition=holy-cow # Partition
#SBATCH -n 8 # Number of Nodes required
#SBATCH -N 1
#SBATCH --mem=96G # Memory request 
#SBATCH -t 14-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o minimap_%J_%a.out # Standard output
#SBATCH -e minimap_%J_%a.err # Standard error

source ~/mambaforge/bin/activate minimap

export SLURM_TMPDIR=/n/holyscratch01/informatics/dkhost/tmp
export TMPDIR=/n/holyscratch01/informatics/dkhost/tmp

#Assign assembly and read files from spec sheet
SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" samples.txt)

ID=$(echo "$SAMPLE" | cut -f1)
READS=$(echo "$SAMPLE" | cut -f2)

REF="/n/holylfs05/LABS/informatics/Lab/scrubjay/GENOMIC_RESOURCES/assemblies/reference/results/FINAL/aphWoo.v1.fa"

#minimap2 -L -t 8 -ax map-hifi $REF $READS > ${ID}_vs_aphWoo1Ref.sam

source ~/mambaforge/bin/activate samtools
samtools view -T $REF -h -b -F 4 ${ID}_vs_aphWoo1Ref.sam | samtools sort -@ 8 -o ${ID}_vs_aphWoo1Ref.mapped.sorted.bam - 
samtools index ${ID}_vs_aphWoo1Ref.mapped.sorted.bam
