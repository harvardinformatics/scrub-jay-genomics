#!/bin/bash
#SBATCH -J map_sats # A single job name for the array
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


SAMPLE="SJ_49358"
REF="/n/holylfs04/LABS/edwards_lab/Lab/sedwards/scrubjay/assemblies/new_AC_birds_fastqs/SJ_49358_cat_all_new.fasta.gz"
SATS="/n/holylfs05/LABS/informatics/Lab/scrubjay/GENOMIC_RESOURCES/repeat_annotation/results/sj_sats_combined_assem.fasta"
mkdir sj_sats_vs_${SAMPLE}/
cd sj_sats_vs_${SAMPLE}/

####COMMANDS####

source ~/mambaforge/bin/activate minimap
#Map vs GENOME
minimap2 -t 16 -c -N1000000 -f1000 -r100,100 <($SRFUTILS enlong ${SATS}) $REF > sj_sats_srf-aln_vs_${SAMPLE}.paf
$SRFUTILS paf2bed sj_sats_srf-aln_vs_${SAMPLE}.paf > sj_sats_srf-aln_vs_${SAMPLE}.bed
$SRFUTILS bed2abun sj_sats_srf-aln_vs_${SAMPLE}.bed > sj_sats_srf-aln_vs_${SAMPLE}.len

cd ../

###############################
SAMPLE="SJ_49710"
REF="/n/holylfs04/LABS/edwards_lab/Lab/sedwards/scrubjay/assemblies/new_AC_birds_fastqs/SJ_49710_cat_all_new.fasta.gz"
SATS="/n/holylfs05/LABS/informatics/Lab/scrubjay/GENOMIC_RESOURCES/repeat_annotation/results/sj_sats_combined_assem.fasta"
mkdir sj_sats_vs_${SAMPLE}/
cd sj_sats_vs_${SAMPLE}/

####COMMANDS####

source ~/mambaforge/bin/activate minimap
#Map vs GENOME
minimap2 -t 16 -c -N1000000 -f1000 -r100,100 <($SRFUTILS enlong ${SATS}) $REF > sj_sats_srf-aln_vs_${SAMPLE}.paf
$SRFUTILS paf2bed sj_sats_srf-aln_vs_${SAMPLE}.paf > sj_sats_srf-aln_vs_${SAMPLE}.bed
$SRFUTILS bed2abun sj_sats_srf-aln_vs_${SAMPLE}.bed > sj_sats_srf-aln_vs_${SAMPLE}.len
