
#######################################
# Structural Variant Analysis with SyRI
#######################################

# This section performs structural variant analysis using SyRI.
# It includes scaffolding with RagTag, mapping with Minimap2, and variant calling with SyRI.

#####################################
# Scaffolding Assemblies using RagTag
#####################################

# Generate RagTag scaffolding commands for each sample and append to a script
for a in "${samples[@]}"; do
    echo "
########
# ${a}
########
ragtag.py scaffold \$reference_genome \$FASTA_${a}_hap1 -o ${a}_hap1_ragtag -t 20
ragtag.py scaffold \$reference_genome \$FASTA_${a}_hap2 -o ${a}_hap2_ragtag -t 20
" >> Codes_RagTag_haps.sh
done

########################################################
# Mapping RagTag Assemblies to VGP Genome using Minimap2
########################################################

# Generate minimap2 mapping commands for SyRI and append to a script
for a in "${samples[@]}"; do
    echo "
########
# ${a}
########
minimap2 -ax asm5 reference_genome_autosomes.fa ${a}_hap1_ragtag.fa -t 34 --eqx | samtools sort -m4G -@ 34 -O BAM -o \$output_folder/${a}_hap1_syri.bam
minimap2 -ax asm5 reference_genome_autosomes.fa ${a}_hap2_ragtag.fa -t 34 --eqx | samtools sort -m4G -@ 34 -O BAM -o \$output_folder/${a}_hap2_syri.bam
" >> Codes_minimap.sh
done

########################################
# Structural Variant Detection with SyRI
########################################

# Generate SyRI commands for each sample and append to a script
for a in "${samples[@]}"; do
    echo "
########
# ${a}
########
samtools index \$output_folder/${a}_hap1_syri.bam
samtools index \$output_folder/${a}_hap2_syri.bam

syri -c \$output_folder/${a}_hap1_syri.bam -r VGP_autosomes_for_SyRi.fa -q ${a}_hap1_autosomes_for_SyRi.fa -F B --prefix VGP_${a}_hap1_ --nc 30 --dir ..
syri -c \$output_folder/${a}_hap2_syri.bam -r VGP_autosomes_for_SyRi.fa -q ${a}_hap2_autosomes_for_SyRi.fa -F B --prefix VGP_${a}_hap2_ --nc 30 --dir ..
" >> Codes_SyRi.sh
done
