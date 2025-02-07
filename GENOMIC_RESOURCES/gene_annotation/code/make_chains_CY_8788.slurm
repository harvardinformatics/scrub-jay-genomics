#!/bin/bash
#SBATCH -n 8                # Number of cores
#SBATCH -N 1
#SBATCH -t 7-0:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p holy-cow   # Partition to submit to
#SBATCH --mem=64G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o LASTZ_CY_8788.out # File to which STDOUT will be written
#SBATCH -e LASTZ_CY_8788.err # File to which STDERR will be written
#SBATCH -J makechains

source ~/mambaforge/bin/activate nextflow

TARID="bGalGal1"
TARGET="/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/REORG/GENOMIC_RESOURCES/gene_annotation/results/chicken_TOGA_output/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic_labelFix.fna"
QUEID="CY_8788"
QUERY="/n/holylfs05/LABS/informatics/Lab/scrubjay/GENOMIC_RESOURCES/gene_annotation/results/suppBirds_TOGA/CY_8788.p_ctg.fa.masked"
OUTDIR="/n/holylfs05/LABS/informatics/Lab/scrubjay/GENOMIC_RESOURCES/gene_annotation/results/suppBirds_TOGA/CY_8788_LASTZ/"


/n/holylfs05/LABS/informatics/Users/dkhost/make_lastz_chains/make_chains.py $TARID $QUEID $TARGET $QUERY --project_dir $OUTDIR --executor slurm --executor_partition holy-cow
