#!/bin/bash
#SBATCH -J srf_test # A single job name for the array
#SBATCH --partition=holy-smokes # Partition
#SBATCH -n 16 # Number of Nodes required
#SBATCH --mem=64G # Memory request 
#SBATCH -t 7-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o test.out # Standard output
#SBATCH -e test.err # Standard error

KMC="/n/holylfs05/LABS/informatics/Users/dkhost/kmc/kmc"
KMC_DUMP="/n/holylfs05/LABS/informatics/Users/dkhost/kmc/kmc_dump"
SRF="/n/holylfs05/LABS/informatics/Users/dkhost/srf/srf"
SRFUTILS="/n/holylfs05/LABS/informatics/Users/dkhost/srf/srfutils.js"

REF="/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/results//HIC_REFERENCE_ASSEMS/AW_365336_FSJragtag.v1.fasta"

#$KMC -fq -k151 -t16 -ci100 -cs100000 @fofn_test.txt count.kmc tmp_dir
#$KMC_DUMP count.kmc count.txt

#$SRF -p test count.txt > srf.fa

source ~/anaconda3/bin/activate minimap
minimap2 -c -N1000000 -f1000 -r100,100 <($SRFUTILS enlong srf.fa) $REF > srf-aln.paf
$SRFUTILS paf2bed srf-aln.paf > srf-aln.bed
$SRFUTILS bed2abun srf-aln.bed > srf-aln.len
