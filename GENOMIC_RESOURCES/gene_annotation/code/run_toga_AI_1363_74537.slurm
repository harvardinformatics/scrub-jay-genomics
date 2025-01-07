#!/bin/bash

#SBATCH -t 14-00:00:00
#SBATCH --mem 16G
#SBATCH -n 1
#SBATCH --partition holy-cow,unrestricted
#SBATCH -J toga

source ~/mambaforge/bin/activate nextflow

CHAIN="/n/holylfs05/LABS/informatics/Lab/scrubjay/GENOMIC_RESOURCES/gene_annotation//results/suppBirds_TOGA/AI_1363_74537_LASTZ/bGalGal1.AI_1363_74537.allfilled.chain.gz"
BED="/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/REORG/GENOMIC_RESOURCES/gene_annotation/results/chicken_TOGA_output/GCF_016699485.2_genomic.FILTERED.bed"
REF2BIT="/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/REORG/GENOMIC_RESOURCES/gene_annotation/results/chicken_TOGA_output/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic_labelFix.2bit"
TAR2BIT="/n/holylfs05/LABS/informatics/Lab/scrubjay/GENOMIC_RESOURCES/gene_annotation//results/suppBirds_TOGA/AI_1363_74537_LASTZ/AI_1363_74537.2bit"
ISO="/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/REORG/GENOMIC_RESOURCES/gene_annotation/results/chicken_TOGA_output/isoforms_GCF_016699485.2_genomic.tsv"

/n/holylfs05/LABS/informatics/Users/dkhost/TOGA/toga.py $CHAIN $BED $REF2BIT $TAR2BIT -i $ISO --nc /n/holylfs05/LABS/informatics/Users/dkhost/TOGA/nextflow_config_files/ --ms --project_name gGal_Toga_AI_1363_74537 --min_score 16000
