We projected the PGGB graph into a VCF on the reference coordinates for analysis using vg and bcftools. The code for this step is in: `export_vcf` subdirectory. We first ran vg deconstruct, then vg wave, then sorted and indexed the resulting vcfs and merged with the `vcf_fix_and_merge_pggb_communities.sh` script.

We also generated PAF files of each haplotype assembly mapped to the reference using minimap, see `minimap` subdirectory.

We also used SyRI to anaylze inversions: `SyRi_inversions.sh`

## VCF Generation

In order to produce a final vcf, we cleaned and normalized the inputs, extracted the CY genotype as the ancestral allele, remerged, and then computed a variant table. Specifically:

First, we fill in missing individuals by adding sample columns with only missing data for each individual that is not present in the input community vcf, and concatenate each individual community vcf together: `vcf_fix_and_merge_pggb_communities.sh`

Second, we get the ancestral allele:
bcftools view -s CY_8788 -Oz -o pggb_cleaned_outgrouponly.vcf.gz --write-index pggb_cleaned_final.vcf.gz
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' pggb_cleaned_outgrouponly.vcf.gz > aa.tab
python calc_anc_allele.py
bgzip aa.processed.tab 
tabix -s1 -b2 -e2 aa.processed.tab.gz 
echo `##INFO=<ID=AA,Number=1,Type=Character,Description="Ancestral allele">` > hdr.txt
bcftools annotate -x INFO/CONFLICT,INFO/LV,INFO/PS,INFO/AT,INFO/TYPE,INFO/LEN,INFO/ORIGIN,INFO/NS pggb_cleaned_final.vcf.gz | bcftools annotate -a aa.processed.tab.gz -c CHROM,POS,REF,ALT,.INFO/AA -h hdr.txt -Oz -o pggb_recode_aa.vcf.gz --write-index

This sets INFO/AA on all records, with missing ancestral information as "."

The, we do some final processing:
bcftools view -Ou -a -s ^CY_8788 pggb_recode_aa.vcf.gz | bcftools annotate -Ou -x INFO/AF,INFO/F_MISSING | bcftools view -c 1 -Ou | bcftools annotate -a genomic_overlaps.bed -h gr_header.txt -c CHROM,FROM,TO,GR --merge-logic GR:unique -Oz -o pggb_ingroup_only.vcf.gz --write-index
bcftools query -f "%CHROM\t%POS0\t%END\t%TYPE\t%INFO/GR\t%REF\t%ALT\t%INFO/AA\t%INV[\t%GT]" pggb_ingroup_only.vcf.gz > pggb_variation.tab 
cut -f1,2,3 pggb_variation.tab | bedtools intersect -a - -b ../../GENOMIC_RESOURCES/genomic_overlaps/genomic_overlaps.bed -loj | cut -f 1,2,3,7 > pggb_variation_genomic_overlaps.tab
cut -f1,2,3 pggb_variation.tab | bedtools intersect -a - -b ../../GENOMIC_RESOURCES/genomic_overlaps/aphWoo_repeats.bed -loj | cut -f 1,2,3,7 > pggb_variation_repeat_overlaps.tab

Finally, we process the tab file into an analysis-ready output with `compute_var_table.py`

(these details are included in `fix_merged_pggb_vcf.sh` and `process_vcfs.sh`)

