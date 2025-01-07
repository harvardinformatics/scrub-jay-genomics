#!/bin/bash
#SBATCH -J repMask # A single job name for the array
#SBATCH -a 1-3
#SBATCH --partition=holy-cow,intermediate # Partition
#SBATCH -n 8 # Number of Nodes required
#SBATCH --mem=64G # Memory request 
#SBATCH -t 7-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o repmask_%J_%a.out # Standard output
#SBATCH -e repmask_%J_%a.err # Standard error

REPLIB="/n/holylfs05/LABS/informatics/Lab/scrubjay/GENOMIC_RESOURCES/repeat_annotation/results/AW_365336_combined_repeats_v2.fasta"
FASTA=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /n/holylfs05/LABS/informatics/Lab/scrubjay/GENOMIC_RESOURCES/repeat_annotation/code/suppBirds_fofn.txt)
SAMPLE=$(basename $FASTA .fa)

mkdir -p /n/holylfs05/LABS/informatics/Lab/scrubjay/GENOMIC_RESOURCES/repeat_annotation/results/RepeatMasker/${SAMPLE}
cd /n/holylfs05/LABS/informatics/Lab/scrubjay/GENOMIC_RESOURCES/repeat_annotation/results/RepeatMasker/${SAMPLE}

/n/holyscratch01/informatics/dkhost/holylfs_shortcut/dfam-tetools.sh --container /n/holyscratch01/informatics/dkhost/holylfs_shortcut/dfam-tetools-latest.sif -- RepeatMasker -xsmall -pa 8 -gff -no_is -dir . -lib $REPLIB $FASTA
