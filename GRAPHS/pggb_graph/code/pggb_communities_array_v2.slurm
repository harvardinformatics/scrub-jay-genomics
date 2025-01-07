#!/bin/bash
#SBATCH -J pggb # A single job name for the array
#SBATCH -a 1
#SBATCH --partition=holy-cow # Partition
#SBATCH -n 12 # Number of Nodes required
#SBATCH --mem=250G # Memory request 
#SBATCH -t 45-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o pggb_%A_%a.out # Standard output
#SBATCH -e pggb_%A_%a.err # Standard error

#source ~/anaconda3/bin/activate pggb
#module load samtools
DIR="/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/results/PGGB"

COM=$(sed "${SLURM_ARRAY_TASK_ID}q;d" results/PGGB/chrom_communities_redo_fofn.txt)
SAMPLE=$(basename $COM .fa)

singularity exec pggb_native/pggb_cannon.sif pggb -i $COM -o ${DIR}/${SAMPLE} -p 94 -s 100000 -n 91 -t 12
