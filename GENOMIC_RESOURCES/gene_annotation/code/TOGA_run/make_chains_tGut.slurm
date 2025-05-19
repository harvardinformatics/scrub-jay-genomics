#!/bin/bash
#SBATCH -n 8                # Number of cores
#SBATCH -N 1
#SBATCH -t 7-0:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared,holy-cow   # Partition to submit to
#SBATCH --mem=120G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o outfile2.out # File to which STDOUT will be written
#SBATCH -e outfile2.err # File to which STDERR will be written
#SBATCH -J makechains

source ~/mambaforge/bin/activate nextflow

TARID="tGut"
TARGET="GCF_003957565.2_bTaeGut1.4.pri_genomic.fna"
QUEID="AW_365336"
QUERY="AW_365336_FSJragtag.v1.ALTLABEL.fasta.masked"
OUTDIR="/n/holyscratch01/informatics/dkhost/bird_TOGA/scrubjay_vs_tgut/chains"


#/n/holylfs05/LABS/informatics/Users/dkhost/make_lastz_chains/make_chains.py $TARID $QUEID $TARGET $QUERY --project_dir $OUTDIR --executor slurm --executor_partition shared --seq1_chunk 50000000 --seq2_chunk 10000000
/n/holylfs05/LABS/informatics/Users/dkhost/make_lastz_chains/make_chains.py $TARID $QUEID $TARGET $QUERY --project_dir $OUTDIR --continue_arg cat --executor local --seq1_chunk 50000000 --seq2_chunk 10000000
