#!/bin/bash

#fix ploidy first
bcftools +fixploidy pggb_merged.vcf.gz > pggb_merged_fixploidy.vcf

#fill tags
bcftools +fill-tags pggb_merged_fixploidy.vcf -Ob -o pggb_merged_cleanup.bcf -- -t AN,AC,AF

#normalized
#this step seems to introduce some INV=1 INFO lines in addition to INV=YES
bcftools norm -m+both -f ../../GENOMIC_RESOURCES/assemblies/reference/results/FINAL/aphWoo.v1.fa pggb_merged_cleanup.bcf > pggb_merged_normalized.vcf

#clean
bcftools sort pggb_merged_normalized.vcf | bcftools +fill-tags -- -t AN,AC,AF,F_MISSING |  bcftools view -o pggb_cleaned_final.vcf.gz -Oz

#biallelic snps
 bcftools view -m2 -M2 -v snps -o pggb_cleaned_final_biallelic_snps.vcf.gz -Oz --write-index pggb_cleaned_final.vcf.gz 
