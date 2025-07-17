#!/bin/bash
#SBATCH -J ragtag_array # A single job name for the array
#SBATCH -a 1-3
#SBATCH --partition=holy-smokes # Partition
#SBATCH -n 8 # Number of Nodes required
#SBATCH --mem=48G # Memory request 
#SBATCH -t 7-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o ragtag_%A_%a.out # Standard output
#SBATCH -e ragtag_%A_%a.err # Standard error

source ~/anaconda3/bin/activate ragtag

#Assign assembly and read files from spec sheet
INFO=$(sed "${SLURM_ARRAY_TASK_ID}q;d" hic_fofn2.txt)

SAMPLE=$(echo $INFO | awk 'BEGIN{FS=" "} {print $1}')

ASSEM="/n/holyscratch01/informatics/dkhost/scrub-jay-genomics_shortcut/workflow/results/HiC_scaffolding_YaHS/${SAMPLE}/${SAMPLE}_yahs_scaff_scaffolds_final.fa"

REF="/n/holylfs05/LABS/informatics/Everyone/scrubjay/RAW_DATA/bCorHaw1.pri.cur.20210920.fasta"

OUT="/n/holyscratch01/informatics/dkhost/scrub-jay-genomics_shortcut/workflow/results/ragtag_bCorHaw"

#minimap2 --cs -cx asm20 -t 8 $REF $ASSEM > ${OUT}/${SAMPLE}_vs_bCorHaw.paf 

ragtag.py scaffold -o ${OUT}/${SAMPLE}_ragtag --mm2-params='-t 8 -x asm10' $REF $ASSEM
