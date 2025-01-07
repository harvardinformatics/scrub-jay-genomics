#!/bin/bash

#SBATCH -t 14-00:00:00
#SBATCH --mem 16G
#SBATCH -n 1
#SBATCH --partition holy-cow,unrestricted
#SBATCH -J toga

source ~/mambaforge/bin/activate nextflow

CHAIN="/n/holyscratch01/informatics/dkhost/bird_TOGA/scrubjay_vs_tgut/chains/tGut.AW_365336.allfilled.chain.gz"
BED="/n/holyscratch01/informatics/dkhost/bird_TOGA/scrubjay_vs_tgut/GCF_003957565.2_bTaeGut1.4.pri_genomic.FILTERED.bed"
REF2BIT="/n/holyscratch01/informatics/dkhost/bird_TOGA/scrubjay_vs_tgut/GCF_003957565.2_bTaeGut1.4.pri_genomic.2bit"
TAR2BIT="/n/holyscratch01/informatics/dkhost/bird_TOGA/scrubjay_vs_tgut/AW_365336_FSJragtag.v1.ALTLABEL.fasta.masked.2bit"
ISO="/n/holyscratch01/informatics/dkhost/bird_TOGA/scrubjay_vs_tgut/isoforms_GCF_003957565.2_bTaeGut1.4.pri_genomic.tsv"

/n/holylfs05/LABS/informatics/Users/dkhost/TOGA/toga.py $CHAIN $BED $REF2BIT $TAR2BIT -i $ISO --nc /n/holylfs05/LABS/informatics/Users/dkhost/TOGA/nextflow_config_files/ --ms
