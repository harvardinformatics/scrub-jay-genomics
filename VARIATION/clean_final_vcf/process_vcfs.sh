bcftools view -s CY_8788 -Oz -o pggb_cleaned_outgrouponly.vcf.gz --write-index pggb_cleaned_final.vcf.gz
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' pggb_cleaned_outgrouponly.vcf.gz > aa.tab
python calc_anc_allele.py
bgzip aa.processed.tab 
tabix -s1 -b2 -e2 aa.processed.tab.gz 
echo '##INFO=<ID=AA,Number=1,Type=Character,Description="Ancestral allele">' > hdr.txt
bcftools annotate -x INFO/CONFLICT,INFO/LV,INFO/PS,INFO/AT,INFO/TYPE,INFO/LEN,INFO/ORIGIN,INFO/NS pggb_cleaned_final.vcf.gz | bcftools annotate -a aa.processed.tab.gz -c CHROM,POS,REF,ALT,.INFO/AA -h hdr.txt -Oz -o pggb_recode_aa.vcf.gz --write-index

## this sets INFO/AA on all records, with missing ancestral information as "."

Final processing
bcftools view -Ou -a -s ^CY_8788 pggb_recode_aa.vcf.gz | bcftools annotate -Ou -x INFO/AF,INFO/F_MISSING | bcftools view -c 1 -Ou | bcftools annotate -a genomic_overlaps.bed -h gr_header.txt -c CHROM,FROM,TO,GR --merge-logic GR:unique -Oz -o pggb_ingroup_only.vcf.gz --write-index
bcftools query -f "%CHROM\t%POS0\t%END\t%TYPE\t%INFO/GR\t%REF\t%ALT\t%INFO/AA\t%INV[\t%GT]" pggb_ingroup_only.vcf.gz > pggb_variation.tab 

# generate a check on genomic overlaps
cut -f1,2,3 pggb_variation.tab | bedtools intersect -a - -b ../../GENOMIC_RESOURCES/genomic_overlaps/genomic_overlaps.bed -loj | cut -f 1,2,3,7 > pggb_variation_genomic_overlaps.tab
cut -f1,2,3 pggb_variation.tab | bedtools intersect -a - -b ../../GENOMIC_RESOURCES/genomic_overlaps/aphWoo_repeats.bed -loj | cut -f 1,2,3,7 > pggb_variation_repeat_overlaps.tab
