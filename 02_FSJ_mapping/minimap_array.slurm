#!/bin/bash
#SBATCH -J minimap_array # A single job name for the array
#SBATCH -a 1-11
#SBATCH --partition=holy-smokes # Partition
#SBATCH -n 8 # Number of Nodes required
#SBATCH --mem=48G # Memory request 
#SBATCH -t 7-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o minimap_%A_%a.out # Standard output
#SBATCH -e minimap_%A_%a.err # Standard error

source activate minimap

p=$(sed "${SLURM_ARRAY_TASK_ID}q;d" assemblies_fofn.txt)
sample=$(basename $p .fa)

ref="/n/holylfs05/LABS/informatics/Everyone/scrubjay/RAW_DATA/FSJ.v2_1.fasta"

out="/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/results/FSJ_mapping/"

minimap2 --cs -cx asm20 -t 8 $ref $p > ${out}/${sample}_vs_FSJ.paf 

sort -k6,6 -k8,8n ${out}/${sample}_vs_FSJ.paf | paftools.js call -f $ref -L10000 -l1000 - > ${out}/${sample}_vs_FSJ.vcf
