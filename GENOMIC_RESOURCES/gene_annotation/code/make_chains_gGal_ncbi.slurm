#!/bin/bash
#SBATCH -n 8                # Number of cores
#SBATCH -N 1
#SBATCH -t 7-0:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p holy-cow   # Partition to submit to
#SBATCH --mem=64G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o makechains2.out # File to which STDOUT will be written
#SBATCH -e makechains2.err # File to which STDERR will be written
#SBATCH -J makechains

source ~/mambaforge/bin/activate nextflow

TARID="bGalGal1"
TARGET="GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic_labelFix.fna"
QUEID="AW_365336"
QUERY="AW_365336_FSJragtag.v1.ALTLABEL.fasta.masked"
OUTDIR="/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/REORG/GENOMIC_RESOURCES/gene_annotation/results/other/chains/"


#/n/holylfs05/LABS/informatics/Users/dkhost/make_lastz_chains/make_chains.py $TARID $QUEID $TARGET $QUERY --project_dir $OUTDIR --executor slurm --executor_partition shared

/n/holylfs05/LABS/informatics/Users/dkhost/make_lastz_chains/make_chains.py $TARID $QUEID $TARGET $QUERY --project_dir $OUTDIR --continue_arg cat --executor local
