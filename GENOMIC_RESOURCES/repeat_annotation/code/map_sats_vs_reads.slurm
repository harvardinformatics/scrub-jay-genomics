#!/bin/bash
#SBATCH -J srf_test # A single job name for the array
#SBATCH -a 1-45
#SBATCH --partition=holy-smokes # Partition
#SBATCH -n 8 # Number of Nodes required
#SBATCH --mem=64G # Memory request 
#SBATCH -t 3-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o sats_vs_reads_%J.out # Standard output
#SBATCH -e sats_vs_reads_%J.err # Standard error

KMC="/n/holylfs05/LABS/informatics/Users/dkhost/kmc/kmc"
KMC_DUMP="/n/holylfs05/LABS/informatics/Users/dkhost/kmc/kmc_dump"
SRF="/n/holylfs05/LABS/informatics/Users/dkhost/srf/srf"
SRFUTILS="/n/holylfs05/LABS/informatics/Users/dkhost/srf/srfutils.js"


p=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/sample_sheet.txt)
SAMPLE=$(echo -n "$p" | cut -f1)

SATS="/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/results/satellite/sj_sats_combined_assem.fa"

TEMP1=$(echo -n "$p" | cut -f3 | tr ',' ' ')

echo $SAMPLE
echo $TEMP1

#mkdir /n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/results/satellite/${SAMPLE}_vs_reads/
cd /n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/results/satellite/sj_sats_vs_${SAMPLE}/


####COMMANDS####
#zcat $TEMP1 > tmp_combined_reads.fastq

source ~/anaconda3/bin/activate minimap
#Map vs READS
#minimap2 -t 8 -c -N1000000 -f1000 -r100,100 <($SRFUTILS enlong ${SATS}) tmp_combined_reads.fastq > sj_sats_vs_${SAMPLE}_reads.paf
$SRFUTILS paf2bed sj_sats_vs_${SAMPLE}_reads.paf > sj_sats_vs_${SAMPLE}_reads.bed
$SRFUTILS bed2abun sj_sats_vs_${SAMPLE}_reads.bed > sj_sats_vs_${SAMPLE}_reads.len

rm -f tmp_combined_reads.fastq
