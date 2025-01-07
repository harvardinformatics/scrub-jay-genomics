#!/bin/bash
#SBATCH -n 12                # Number of cores
#SBATCH -N 1
#SBATCH -t 3-0:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared   # Partition to submit to
#SBATCH --mem=64G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o repmask.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e repmask.err  # File to which STDERR will be written, %j inserts jobid

cd /n/holylfs05/LABS/informatics/Lab/scrubjay/GENOMIC_RESOURCES/repeat_annotation/results/SV_repeatAnnot/

/n/holyscratch01/informatics/dkhost/holylfs_shortcut/dfam-tetools.sh --container /n/holyscratch01/informatics/dkhost/holylfs_shortcut/dfam-tetools-latest.sif -- RepeatMasker -xsmall -pa 9 -gff -no_is -lib /n/holylfs05/LABS/informatics/Lab/scrubjay/GENOMIC_RESOURCES/repeat_annotation/results/AW_365336_combined_repeats_v2.fasta /n/holylfs05/LABS/informatics/Lab/scrubjay/GENOMIC_RESOURCES/repeat_annotation/results/SV_repeatAnnot/pggb_merged_INDEL_alleles.short.fasta
