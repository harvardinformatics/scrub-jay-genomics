#!/bin/bash
#SBATCH -J busco_array # A single job name for the array
#SBATCH -a 1-11
#SBATCH --partition=holy-smokes # Partition
#SBATCH -n 4 # Number of Nodes required
#SBATCH --mem=24G # Memory request 
#SBATCH -t 1-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o busco_%A_%a.out # Standard output
#SBATCH -e busco_%A_%a.err # Standard error

source activate busco

p=$(sed "${SLURM_ARRAY_TASK_ID}q;d" assemblies_fofn.txt)

sample=$(basename $p .fa.gz)

#gunzip $p > ${sample}.fa

busco -f -c 4 -m geno -i $p -o ${sample}_busco -l /n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/01_assembly_qc/aves_odb10/

#rm ${sample}.fa
