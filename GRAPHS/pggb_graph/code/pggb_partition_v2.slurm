#!/bin/bash
#SBATCH -J pggb_part # A single job name for the array
#SBATCH --partition=holy-cow # Partition
#SBATCH -n 16 # Number of Nodes required
#SBATCH --mem=400G # Memory request 
#SBATCH -t 45-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o pggb_%J.out # Standard output
#SBATCH -e pggb_%J.err # Standard error

source ~/anaconda3/bin/activate pggb
module load samtools
DIR="/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/results/PGGB"

#samtools faidx ${DIR}/combined_assemblies.fa
singularity exec pggb_1656841325.sif wfmash ${DIR}/combined_assemblies.fasta -p 94 -n 90 -t 16 -m -w 1024 > ${DIR}/combined_assemblies.partition.paf
