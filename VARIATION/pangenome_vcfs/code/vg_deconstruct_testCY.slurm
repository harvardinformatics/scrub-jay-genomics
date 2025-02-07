#!/bin/bash
#SBATCH -J vg_decon # A single job name for the array
#SBATCH --partition=shared # Partition
#SBATCH -n 4 # Number of Nodes required
#SBATCH --mem=32G # Memory request 
#SBATCH -t 3-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o vgde_%A_%a.out # Standard output
#SBATCH -e vgde_%A_%a.err # Standard error


COM="/n/holylfs05/LABS/informatics/Lab/scrubjay/GRAPHS/pggb_graph/results/pggb_graphs/allbird_community.21/allbird_community.21.final.gfa"
mkdir -p testCY

/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/vg deconstruct -P CY_8788 -H '#' -e -a -t 4 $COM  > testCY/community21_testCY.vcf
