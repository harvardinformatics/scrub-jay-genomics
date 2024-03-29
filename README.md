# Scrub Jay Genomics

This is a minimal initial repository to track code for diploid assembly and long-read population genomics of scrub jays. 

## Where are the results?
Revelant output files can be found on the Cannon cluster at `/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/results`

## Results subdirectories
### REFERENCE ASSEMBLIES:  
- `HIC_REFERENCE_ASSEMS/AW_365336_FSJragtag.v1.fasta`: scaffolded (HiC + Ragtag) *A. woodhouseii* reference female  
- `HIC_REFERENCE_ASSEMS/AW_365336_FSJragtag.v1.ALTLABEL.fasta`: as above, but with headers in Panspec format (for PGGB)  
- `HIC_REFERENCE_ASSEMS/AW_365336_combined_repeats_v2.fasta`: repeat library used for annotation (RepeatModeler + SRF)  
- `HIC_REFERENCE_ASSEMS/AW_365336_FSJragtag.v1.fasta.tbl`: tabular output of RepeatMasker for reference individual  
- `HIC_REFERENCE_ASSEMS/AW_365336_FSJragtag.v1.fasta.out.gff`: RepeatMasker annotation  
- `HIC_REFERENCE_ASSEMS/AW_365338_FSJragtag.v1.fasta`: scaffolded reference for *A. woodhouseii* reference male   
### INDIVIDUAL ASSEMBLIES  
- `assemblies/*.p_ctg.fasta`: unscaffolded primary assemblies  
    - AC = *A. coerulescens*, Florida scrub jay 
    - AI = *A. insularis*, Island scrub jay  
    - AW = *A. woodhouseii*, Woodhouse scrub jay  
    - CY = Yucatan scrub jay  
- `assemblies/*hap[1|2].p_ctg.fasta`: unscaffolded haplotype assemblies  
- `assembly_qc/`: basic assembly stats  
- `assembly_qc/ASSEMBLY_STATS.tsv`: summary file of basic stats, primary and haplotype  
- `assembly_qc/READS_STATS.tsv`: summary file with basic stats for FASTQ files used for each assembly
### PANGENOME GRAPHS AND ANALYSIS  
- *NB:* as of 2/1/23, several communities are still in process of construction
- `PGGB/combined_assemblies.partition.paf`: alignment file of *all assembly haplotypes*, plus reference and Yucatan jay
- `PGGB/communities/`: lists of all communities partitioned by wfmash  
- `PGGB/allbird_community.[0-9]`: communities containing **reference chromosomes**  
    - `PGGB/allbird_community.[0-9]/*.paf`: alignment file of sequences in community  
    - `PGGB/allbird_community.[0-9]/*.og`: PGGB graph format of alignment (input to ODGI for visualization)  
    - `PGGB/allbird_community.[0-9]/*.gfa`: standard graph format of alignment  
    - `PGGB/allbird_community.[0-9]/*.nameFix.vcf.gz`: varaiant call format file deconstructed from .og file using vg deconstruct. *Note: 'nameFix' version has the reference geneome ID as 'aphWoo1' to properly recolve haplotypes*  
    - `PGGB/allbird_community.[0-9]/*final_nameFix_bub_wave.vcf`: normalized and deconvoluted VCF file (run thru vcfwave and vcfbub). *Use this file for pop gen analysis!*  
    - `PGGB/allbird_community.[0-9]/*bub_wave_A[W|I|C]_bialle.vcf`: normalized and deconvoluted VCF file of only biallelic SNPs, split by species
- `PGGB/allbird_community.[0-9]_unplaced`: communities containing **unplaced reference scaffolds**  
    - Contains same files as above   

### SRF SATELLITES  
- `satellite/sj_sats_combined_assem.fa`: KMC+SRF output with satellites from all combined (primary) assemblies  
- `satellite/sj_sats_vs_*/`: combined satellites mapped against individual assemblies  
    - `sj_sats_srf-aln_vs_A[W|I|C]_*.bed`: combined satellites mapped vs the individual genome  
    - `sj_sats_srf-aln_vs_A[W|I|C]_*_reads.bed`: combined satellites mapped vs the individual genomic **reads**  
        - There are also `.paf` (alignment) files and `.len` (repeat count summary) files for each sample  
- `results/satellite/individual_vs_reads/`: KMC+SRF output from **individual** reads (i.e. NOT combined satellites)  

### GENE ANNOTATION  
- `gene_annotation/stringtie_RNAseq`: annotation with Illumina, done using stringtie + hisat. BUSCO completeness 97%  
    - `scrubjay_stringtie_merge.gtf`: merged GTF file from all samples  
    - `scrubjay_stringtie_transcripts.fa`: transcripts extracted using gffread  
- `gene_annotation/IsoQuant/`: annotation using PacBio Isoseq, done using IsoQuant. BUSCO completeness 88%  
    - `isoquant_merge.gtf`: merged GTF file  
    - `transcripts.fa`: transcripts extracted using gffread