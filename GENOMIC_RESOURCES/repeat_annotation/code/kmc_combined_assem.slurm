#!/bin/bash
#SBATCH -J srf_test # A single job name for the array
#SBATCH --partition=holy-cow # Partition
#SBATCH -n 16 # Number of Nodes required
#SBATCH --mem=400G # Memory request 
#SBATCH -t 7-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o kmc.out # Standard output
#SBATCH -e kmc.err # Standard error

KMC="/n/holylfs05/LABS/informatics/Users/dkhost/kmc/kmc"
KMC_DUMP="/n/holylfs05/LABS/informatics/Users/dkhost/kmc/kmc_dump"
SRF="/n/holylfs05/LABS/informatics/Users/dkhost/srf/srf"
SRFUTILS="/n/holylfs05/LABS/informatics/Users/dkhost/srf/srfutils.js"

REF="/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/results/assemblies/combined_pri_assems.fasta"

$KMC -fm -k151 -t16 -ci100 -cs100000 $REF count.kmc tmp_dir
#echo "KMC DONE"

$KMC_DUMP count.kmc count.txt
echo "KMC DUMP DONE"

$SRF -p sj_sat count.txt > sj_sats_combined_assem.fa
echo "SRF DONE"

#source ~/anaconda3/bin/activate minimap
#minimap2 -c -N1000000 -f1000 -r100,100 <($SRFUTILS enlong srf.fa) $REF > srf-aln.paf
#$SRFUTILS paf2bed srf-aln.paf > srf-aln.bed
#$SRFUTILS bed2abun srf-aln.bed > srf-aln.len
