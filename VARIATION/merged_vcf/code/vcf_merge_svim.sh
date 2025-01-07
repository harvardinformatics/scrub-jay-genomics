
while read vcf; do

    basename -s .hap1_diploid_svim.vcf $vcf > temp.txt
    SAMPLE=$(basename $vcf)
    bcftools reheader -s temp.txt $vcf | bcftools sort -O z > merged_vcf/fixed_$SAMPLE.gz
    bcftools index merged_vcf/fixed_$SAMPLE.gz 

done < svim_vcf_fofn.txt

ls merged_vcf/fixed*.vcf.gz > merged_vcf/merge_fofn.txt

bcftools merge -l merged_vcf/merge_fofn.txt -O z -o merged_vcf/svim_merged.vcf.gz