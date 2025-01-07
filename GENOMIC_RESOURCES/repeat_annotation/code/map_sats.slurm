#!/bin/bash
#SBATCH -J map_sats # A single job name for the array
#SBATCH -a 1-45
#SBATCH --partition=holy-smokes # Partition
#SBATCH -n 16 # Number of Nodes required
#SBATCH --mem=64G # Memory request 
#SBATCH -t 7-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o map_sats.out # Standard output
#SBATCH -e map_sats.err # Standard error

KMC="/n/holylfs05/LABS/informatics/Users/dkhost/kmc/kmc"
KMC_DUMP="/n/holylfs05/LABS/informatics/Users/dkhost/kmc/kmc_dump"
SRF="/n/holylfs05/LABS/informatics/Users/dkhost/srf/srf"
SRFUTILS="/n/holylfs05/LABS/informatics/Users/dkhost/srf/srfutils.js"


p=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/sample_sheet.txt)
SAMPLE=$(echo -n "$p" | cut -f1)

REF="/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/results/assemblies/${SAMPLE}.p_ctg.fa"
SATS="/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/results/satellite/sj_sats_combined_assem.fa"

mkdir /n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/results/satellite/sj_sats_vs_${SAMPLE}/
cd /n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/results/satellite/sj_sats_vs_${SAMPLE}/

####COMMANDS####

source ~/anaconda3/bin/activate minimap
#Map vs GENOME
minimap2 -c -N1000000 -f1000 -r100,100 <($SRFUTILS enlong ${SATS}) $REF > sj_sats_srf-aln_vs_${SAMPLE}.paf
$SRFUTILS paf2bed sj_sats_srf-aln_vs_${SAMPLE}.paf > sj_sats_srf-aln_vs_${SAMPLE}.bed
$SRFUTILS bed2abun sj_sats_srf-aln_vs_${SAMPLE}.bed > sj_sats_srf-aln_vs_${SAMPLE}.len

