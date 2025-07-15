#!/bin/bash
#SBATCH -J mapHaplo # A single job name for the array
#SBATCH -a 1-45
#SBATCH --partition=holy-cow # Partition
#SBATCH -n 8 # Number of Nodes required
#SBATCH --mem=64G # Memory request 
#SBATCH -t 7-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o map_%A_%a.out # Standard output
#SBATCH -e map_%A_%a.err # Standard error

source ~/mambaforge/bin/activate minimap

DIR="/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/results/haplos_vs_HiC_Ref"

SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" sample_IDs.txt)

REF="/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/results/HIC_REFERENCE_ASSEMS/AW_365336_FSJragtag.v1.ALTLABEL.fasta"
QUERY1="results/assemblies/${SAMPLE}.hap1.p_ctg.fa"
QUERY2="results/assemblies/${SAMPLE}.hap2.p_ctg.fa"

minimap2 -c -L -t 8 -x asm20 $REF $QUERY1 > $DIR/${SAMPLE}.hap1.p_ctg_vs_ref.paf
minimap2 -c -L -t 8 -x asm20 $REF $QUERY2 > $DIR/${SAMPLE}.hap2.p_ctg_vs_ref.paf
