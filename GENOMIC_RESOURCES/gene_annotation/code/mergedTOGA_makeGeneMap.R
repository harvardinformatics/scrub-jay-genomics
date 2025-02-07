library(tidyverse)

trans <- read_tsv("/n/holylfs05/LABS/informatics/Lab/scrubjay/GENOMIC_RESOURCES/gene_annotation/results/merged_TOGA/WORKFLOW/scrubjay_TOGA_merged_sortedTranscripts.bed",col_names=c("scaff","start","end","transID"))

regGene <- read_tsv("/n/holylfs05/LABS/informatics/Lab/scrubjay/GENOMIC_RESOURCES/gene_annotation/results/merged_TOGA/WORKFLOW/scrubjay_TOGA_merged_sortedGenes.bed",col_names=c("scaff","start","end","regGeneID"))

refGene <- read_tsv("/n/holylfs05/LABS/informatics/Lab/scrubjay/GENOMIC_RESOURCES/gene_annotation/results/merged_TOGA/WORKFLOW/scrubjay_TOGA_merged_sortedRefGenes.bed",col_names=c("scaff","start","end","refGeneID"))


mergedGenes <- full_join(refGene, regGene, join_by(scaff, start, end))


by <- join_by(scaff, closest(start >= start), closest(end <= end))
geneMap <- full_join(trans, mergedGenes, by) 

dups <- geneMap[duplicated(geneMap$transID),]

colnames(geneMap) <- c("scaff","transcript_start", "transcript_end", "transcript_ID","locus_start", "locus_end", "ref_GeneID", "toga_GeneID")

write_tsv(geneMap, "/n/holylfs05/LABS/informatics/Lab/scrubjay/GENOMIC_RESOURCES/gene_annotation/results/merged_TOGA/WORKFLOW/gene_map.tsv")
