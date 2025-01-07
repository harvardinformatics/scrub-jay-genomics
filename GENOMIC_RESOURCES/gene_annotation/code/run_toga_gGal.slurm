#!/bin/bash

#SBATCH -t 14-00:00:00
#SBATCH --mem 16G
#SBATCH -n 1
#SBATCH --partition holy-cow,unrestricted
#SBATCH -J toga

source ~/mambaforge/bin/activate nextflow

CHAIN="/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/REORG/GENOMIC_RESOURCES/gene_annotation/results/other/chains/bGalGal1.AW_365336.allfilled.chain.gz"
BED="/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/REORG/GENOMIC_RESOURCES/gene_annotation/results/other/GCF_016699485.2_genomic.FILTERED.bed"
REF2BIT="/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/REORG/GENOMIC_RESOURCES/gene_annotation/results/other/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic_labelFix.2bit"
TAR2BIT="/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/REORG/GENOMIC_RESOURCES/gene_annotation/results/other/AW_365336_FSJragtag.v1.ALTLABEL.fasta.masked.2bit"
ISO="/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/REORG/GENOMIC_RESOURCES/gene_annotation/results/other/isoforms_GCF_016699485.2_genomic.tsv"

/n/holylfs05/LABS/informatics/Users/dkhost/TOGA/toga.py $CHAIN $BED $REF2BIT $TAR2BIT -i $ISO --nc /n/holylfs05/LABS/informatics/Users/dkhost/TOGA/nextflow_config_files/ --ms --project_name gGal_Toga_redo --min_score 16000
