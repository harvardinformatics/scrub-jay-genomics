#!/bin/bash
#SBATCH -J paf2chain # A single job name for the array
#SBATCH --partition=holy-smokes # Partition
#SBATCH -n 4 # Number of Nodes required
#SBATCH -N 1
#SBATCH -a 1-90
#SBATCH --mem=16G # Memory request 
#SBATCH -t 0-04:00 # Maximum execution time (D-HH:MM)
#SBATCH -o paf2chain_%J.out # Standard output
#SBATCH -e paf2chain_%J.err # Standard error


source ~/mambaforge/bin/activate rust
PATH=$PATH:/n/home04/dkhost/.cargo/bin

SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" paf_fofn.txt)
ID=$(basename $SAMPLE .paf)

awk 'BEGIN{FS="\t";OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$NF}' $SAMPLE > temp_${ID}.paf

wgatools paf2chain -t 4 -o $ID.chain temp_${ID}.paf
