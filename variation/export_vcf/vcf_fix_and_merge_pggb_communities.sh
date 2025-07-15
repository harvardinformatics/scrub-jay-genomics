
all_samples=$(bcftools query -l ../results/allbird_community.11_bub_wave_sort.vcf.gz | sort -u)

#get correct sample order so can run concat later (samples must be in same order!)
bcftools query -l allbird_community.1_bub_wave_sort.vcf.gz > samples.txt

# Loop over each VCF file
for vcf in *bub_wave_sort.vcf.gz; do
    # Extract sample names from the current VCF
    samples_in_vcf=($(bcftools query -l $vcf))
    
    # Compare the lists and find missing samples
    missing_samples=($(echo ${all_samples[@]} ${samples_in_vcf[@]} | tr ' ' '\n' | sort | uniq -u))

    # If there are missing samples, create dummy VCFs for each and merge
    if [ ${#missing_samples[@]} -ne 0 ]; then

        # Extract header from the original VCF PROPERLY THIS TIME, JUST THE FORMAT LINES -DEK
        bcftools view -h $vcf | grep '^##' > header.txt
        
        # Loop over each missing sample and create dummy VCF
        for sample in "${missing_samples[@]}"; do

            #ADD THE FORMAT LINES SO VCF FILE CAN BE RECOGNIZED -DEK
            cat header.txt > dummy_$sample.vcf

            # Add missing sample column to header NO NEED FOR THIS WITH HOW DOING HEADER
            echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sample" >> dummy_$sample.vcf

           # Extract CHROM and POS of the first variant
            first_variant_chr=$(bcftools query -f '%CHROM\n' $vcf | head -1)
            first_variant_pos=$(bcftools query -f '%POS\n' $vcf | head -1)

            # Formulate the region string
            first_variant_region="$first_variant_chr:$first_variant_pos"

           # Add .|. genotype for the first variant position
            bcftools query -f '%CHROM\t%POS\t.\t%REF\t%ALT\t.\t.\t.\tGT\t.|\.\n' -r $first_variant_region $vcf >> dummy_$sample.vcf
            bgzip dummy_$sample.vcf
            bcftools index dummy_$sample.vcf.gz

        done
        
        #AT field causing problems; removing it and making temporary vcf file for merging
	bcftools annotate -O z -x INFO/AT $vcf > temp.vcf.gz
        bcftools index temp.vcf.gz

        bcftools merge temp.vcf.gz dummy_*.vcf.gz | bcftools view -S samples.txt -O z -o fixed_$vcf 

       # Merge the original VCF with all dummy VCFs
        #bcftools merge $vcf dummy_*.vcf.gz -Oz -o fixed_$vcf

        rm dummy_*.vcf.gz
        rm dummy_*.vcf.gz.csi

        rm temp.vcf.gz
        rm temp.vcf.gz.csi
    fi
done
rm header.txt
