#!/bin/bash
#SBATCH -J map_sats # A single job name for the array
#SBATCH -a 1-45
#SBATCH --partition=holy-smokes # Partition
#SBATCH -n 16 # Number of Nodes required
#SBATCH --mem=64G # Memory request 
#SBATCH -t 7-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o map_sats.out # Standard output
#SBATCH -e map_sats.err # Standard error

SRF="/n/holylfs05/LABS/informatics/Users/dkhost/srf/srf"
SRFUTILS="/n/holylfs05/LABS/informatics/Users/dkhost/srf/srfutils.js"


PAF=$(sed "${SLURM_ARRAY_TASK_ID}q;d" srf_paf_fofn.txt)
SAMPLE=$(basename $PAF .paf)
OUTDIR=$(echo $SAMPLE | sed 's/srf-aln_//g')

cd ../results/SRF_satellites/$OUTDIR

####COMMANDS####

source ~/anaconda3/bin/activate minimap
#Map vs GENOME
$SRFUTILS paf2bed $PAF > ${SAMPLE}_redo.bed
$SRFUTILS bed2abun ${SAMPLE}_redo.bed > ${SAMPLE}_redo.len

