`%!in%` = Negate(`%in%`)
library(strex)
library(dplyr)
library(tidyverse)
install.packages("assertr")
library(assertr)

#read in pangene output
#Supplementary Table S12
pangene<-read.table("sj-94-a2.pg.CNV.Rtable",header = T,sep = "\t",comment.char = "@")

#create table of Z-linked genes to subtract from pangene table
#the files below are not in Dryad but are easily constructed from annotation files on Dryad
gene_map<-read_delim("scrubjay_genes.tsv", delim = "\t",
                     col_names = c("chr", "st", "en", "strand", "id_info")) %>%
  mutate(id_info = str_remove_all(id_info, "\"")) %>%
  separate_wider_delim(id_info, ";", names=c("tcons_id", "xloc_id")) %>%
  mutate(tcons_id = str_remove(tcons_id, "transcript_id ")) %>%
  mutate(xloc_id = str_remove(xloc_id, "gene_id "))

gene_info<-read_delim("merged_TOGA_tgutRef_noDup.tsv", delim="\t", 
                      col_names = c("tcons_id", "xloc_id", "galGal_id", "tGut_id", "galGal_toga", "tGut_toga", "equiv_class", "gene_name"))

zchrom<-gene_map %>% filter(chr == "aphWoo1#REF#ScYP8k310HRSCAF43chZ_RagTag")
zgenes<-unique(zchrom$xloc_id)
zlist<-gene_info %>% filter(xloc_id %in% zgenes) %>% select(gene_name) %>% distinct()
zlist2<-str_before_first(str_after_first(zlist$gene_name,"gene-"),"_tGut")
zlist3<-str_before_first(str_after_first(zlist$gene_name[zna],"gene-"),"_gGal")
totz<-c(zlist2,zlist3)
noz<-!is.na(totz)
totz2<-totz[noz]
totz2<-unique(totz2)
totz2[682:684]<-c("reg_9383_gGal","reg_9554_gGal","reg_1559_tGut")
write.table(totz2,"Z_linked_gene_list_684.txt",sep = "\t",quote = F,row.names = F,col.names = F)

#read in table of Z-linked genes
totz<-read.table("Z_linked_gene_list_684.txt",sep = "\t")

#remove Z linked genes from pangene table
pangeneauto<-pangene %>% filter(Gene %!in% totz$V1)

#get total copy number of each gene in each haplotype
pangeneauto$sumAI<-rowSums(pangeneauto[,grep("AI_",colnames(pangeneauto))])
pangeneauto$sumAC<-rowSums(pangeneauto[,grep("AC_",colnames(pangeneauto))])
pangeneauto$sumAW<-rowSums(pangeneauto[,grep("AW_",colnames(pangeneauto))])

#create indexes of whether a gene is reliably polymorphic (> 2 haplotypes varying) in each species
#this filter works only for deletions of genes with '1 0' polymorphisms, which is most genes

AIcols<-grep("AI_1",colnames(pangeneauto))
AWcols<-grep("AW_3",colnames(pangeneauto))
ACcols<-grep("AC_1",colnames(pangeneauto))

for(i in 1:dim(pangeneauto)[1]){
  pangeneauto[i,'polyAI']<-ifelse(max(table(t(pangeneauto[i,AIcols]))) >= 29,FALSE,TRUE)
  pangeneauto[i,'polyAC']<-ifelse(max(table(t(pangeneauto[i,ACcols]))) >= 27,FALSE,TRUE)
  pangeneauto[i,'polyAW']<-ifelse(max(table(t(pangeneauto[i,AWcols]))) >= 29,FALSE,TRUE)
}

#make filter to flag genes with only 1/0 polymorphisms
pred<-in_set(1,0)
i<-2
for(i in 1:dim(pangeneauto)[1]){
x<-as.data.frame(t(pangeneauto[i,ACcols]))
colnames(x)<-"counts"
pangeneauto[i,'AC_0_1']<-ifelse((all(pred(x$counts)) & pangeneauto[i,'num_copy_num_states_AC'] > 1),TRUE,FALSE)
}
for(i in 1:dim(pangeneauto)[1]){
  x<-as.data.frame(t(pangeneauto[i,AIcols]))
  colnames(x)<-"counts"
  pangeneauto[i,'AI_0_1']<-ifelse((all(pred(x$counts)) & pangeneauto[i,'num_copy_num_states_AI'] > 1),TRUE,FALSE)
}
for(i in 1:dim(pangeneauto)[1]){
  x<-as.data.frame(t(pangeneauto[i,AWcols]))
  colnames(x)<-"counts"
  pangeneauto[i,'AW_0_1']<-ifelse((all(pred(x$counts)) & pangeneauto[i,'num_copy_num_states_AW'] > 1),TRUE,FALSE)
}

#add number of states of each gene for each species
for(i in 1:dim(pangeneauto)[1]){
pangeneauto[i,'num_copy_num_states_AI']<-length(unique(unlist(as.vector(pangeneauto[i,AIcols]))))
pangeneauto[i,'num_copy_num_states_AC']<-length(unique(unlist(as.vector(pangeneauto[i,ACcols]))))
pangeneauto[i,'num_copy_num_states_AW']<-length(unique(unlist(as.vector(pangeneauto[i,AWcols]))))
}

#determine number of species in which a gene is polymorphic
#using polyA* criterion (at least two haplotypes different)
for(i in 1:dim(pangeneauto)[1]){
  pangeneauto$num_sp_poly_poly[i]<-length(which(pangeneauto[i,grep("polyA",colnames(pangeneauto))] == TRUE))
}

write.table(pangeneauto,"pangene_auto_table_final_poly_stats.txt",sep = "\t",col.names = T,row.names = F,quote = F)
write.table(pangeneauto[,-c(2:95)],"pangene_auto_no_haps_poly_stats_final_table.txt",sep = "\t",col.names = T,row.names = F,quote = F)

#begin calculation of CNV numbers with different filters
#make table longer, pivoting gene counts per haplotype
pangenepoly3<-pangeneauto %>% pivot_longer(c(2:95),names_to = "haplotype",values_to = "copy_num")


#add species designations
pangenepoly3$species<-str_before_first(pangenepoly3$haplotype,"_")

#add individual designations
pangenepoly3$indiv<-str_before_last_dot(pangenepoly3$haplotype)

#generate summaries of CNV numbers with no filter using num_copy_num_states for each species
#for AW
tab<-pangenepoly3 %>% filter(species == "AW") %>% group_by(Gene,species) %>% distinct(num_copy_num_states_AW) 
table(tab$num_copy_num_states_AW)
sum(table(tab$num_copy_num_states_AW)[2:8],na.rm = T)
#for AC
tab<-pangenepoly3 %>% filter(species == "AC") %>% group_by(Gene,species) %>% distinct(num_copy_num_states_AC) 
table(tab$num_copy_num_states_AC)
sum(table(tab$num_copy_num_states_AC)[2:6])
#for AI
tab<-pangenepoly3 %>% filter(species == "AI") %>% group_by(Gene,species) %>% distinct(num_copy_num_states_AI) 
table(tab$num_copy_num_states_AI)
sum(table(tab$num_copy_num_states_AI)[2:8])
##

#now filter by polyA*
#for AW
tab<-pangenepoly3 %>% filter(species == "AW" & polyAW==TRUE) %>% group_by(Gene,species) %>% distinct(num_copy_num_states_AW) 
table(tab$num_copy_num_states_AW)
sum(table(tab$num_copy_num_states_AW),na.rm = T)
#for AC
tab<-pangenepoly3 %>% filter(species == "AC" & polyAC==TRUE) %>% group_by(Gene,species) %>% distinct(num_copy_num_states_AC) 
table(tab$num_copy_num_states_AC)
sum(table(tab$num_copy_num_states_AC))
#for AI
tab<-pangenepoly3 %>% filter(species == "AI" & polyAI==TRUE) %>% group_by(Gene,species) %>% distinct(num_copy_num_states_AI) 
table(tab$num_copy_num_states_AI)
sum(table(tab$num_copy_num_states_AI))

#now do CNVs in at least 6 haplotypes, using sum across haplotypes as metric
#for AW
tab<-pangenepoly3 %>% filter(species == "AW" & (sumAW <= 24 | sumAW >= 36)) %>% group_by(Gene,species) %>% distinct(num_copy_num_states_AW) 
table(tab$num_copy_num_states_AW)
sum(table(tab$num_copy_num_states_AW))
#for AC
tab<-pangenepoly3 %>% filter(species == "AC" & (sumAC <= 22 | sumAC >= 34)) %>% group_by(Gene,species) %>% distinct(num_copy_num_states_AC) 
table(tab$num_copy_num_states_AC)
sum(table(tab$num_copy_num_states_AC))
#for AI
tab<-pangenepoly3 %>% filter(species == "AI" & (sumAI <= 24 | sumAI >= 36)) %>% group_by(Gene,species) %>% distinct(num_copy_num_states_AI) 
table(tab$num_copy_num_states_AI)
sum(table(tab$num_copy_num_states_AI))

#number of genes with biallelic 0/1 PAV polymorphisms
pangeneauto %>% reframe(sumAI10 = sum(AI_0_1))
pangeneauto %>% reframe(sumAW10 = sum(AW_0_1))
pangeneauto %>% reframe(sumAC10 = sum(AC_0_1))

#calculate above in a different way to check
pangenepoly3 %>% filter(species == "AW" & AW_0_1 == TRUE) %>% ungroup() %>% summarize(count = n_distinct(Gene))
pangenepoly3 %>% filter(species == "AI" & AI_0_1 == TRUE) %>% ungroup() %>% summarize(count = n_distinct(Gene))
pangenepoly3 %>% filter(species == "AC" & AC_0_1 == TRUE) %>% ungroup() %>% summarize(count = n_distinct(Gene))

#above filter with frequency of at least two haplotypes
pangenepoly3 %>% filter(species == "AI" & AI_0_1 == TRUE & polyAI == TRUE) %>% ungroup() %>% summarize(count = n_distinct(Gene))
pangenepoly3 %>% filter(species == "AW" & AW_0_1 == TRUE & polyAW == TRUE) %>% ungroup() %>% summarize(count = n_distinct(Gene))
pangenepoly3 %>% filter(species == "AC" & AC_0_1 == TRUE & polyAC == TRUE) %>% ungroup() %>% summarize(count = n_distinct(Gene))

#above filter with frequency of at least six haplotypes
pangenepoly3 %>% filter(species == "AI" & AI_0_1 == TRUE & sumAI < 25) %>% ungroup() %>% summarize(count = n_distinct(Gene))
pangenepoly3 %>% filter(species == "AW" & AW_0_1 == TRUE & sumAW < 25) %>% ungroup() %>% summarize(count = n_distinct(Gene))
pangenepoly3 %>% filter(species == "AC" & AC_0_1 == TRUE & sumAC < 23) %>% ungroup() %>% summarize(count = n_distinct(Gene))

#count multiplications - definition - more than total number of haplotypes in a species
pangenepoly3 %>% filter(species == "AI" & sumAI > 30) %>% ungroup() %>% summarize(count = n_distinct(Gene))
pangenepoly3 %>% filter(species == "AW" & sumAW > 30) %>% ungroup() %>% summarize(count = n_distinct(Gene))
pangenepoly3 %>% filter(species == "AC" & sumAC > 28) %>% ungroup() %>% summarize(count = n_distinct(Gene))

#above filter with frequency of at least two haplotypes
pangenepoly3 %>% filter(species == "AI" & polyAI == 1 & sumAI > 30) %>% ungroup() %>% summarize(count = n_distinct(Gene))
pangenepoly3 %>% filter(species == "AW" & polyAW == 1 & sumAW > 30) %>% ungroup() %>% summarize(count = n_distinct(Gene))
pangenepoly3 %>% filter(species == "AC" & polyAC == 1 & sumAC > 28) %>% ungroup() %>% summarize(count = n_distinct(Gene))

#above filter with frequency of at least six haplotypes
pangenepoly3 %>% filter(species == "AI" & sumAI > 35) %>% ungroup() %>% summarize(count = n_distinct(Gene))
pangenepoly3 %>% filter(species == "AW" & sumAW > 35) %>% ungroup() %>% summarize(count = n_distinct(Gene))
pangenepoly3 %>% filter(species == "AC" & sumAC > 33) %>% ungroup() %>% summarize(count = n_distinct(Gene))

#now count number of genes with 0/0 genotypes
#add individual identities and copy numbers per individual
pangenepoly3<-pangenepoly3 %>% mutate(indiv = str_before_last_dot(haplotype),sum_indiv = sum(copy_num)) 
pangenepoly3<-pangenepoly3 %>% group_by(Gene,indiv) %>% mutate(copy_num1 = copy_num[1],copy_num2 = copy_num[2],sum_indiv=sum(copy_num1,copy_num2)/2) #%>% distinct()
subsetAphHW<-pangenepoly3 %>% ungroup() %>%
  group_by(Gene,indiv) %>% dplyr::reframe(sum_indiv=sum_indiv,species=species,copy_num1 = copy_num[1],copy_num2 = copy_num[2]) %>% distinct()

#table of genes deleted entirely in each species
totdeltab<-genesumtab %>% filter(species %in% c("AW","AC","AI") & genoBBonly > 0)
write.table(totdeltab,"table_genes_all_00_homozygote_deletions.txt",sep = "\t",col.names = T,row.names = F,quote = F)

#Sum CNV genotypes of all genes regardless of copy number per haplotype
#considers copy numbers >= 1 per haplotype as allele "A"; all others allele B
subsetAphHW<-subsetAphHW %>% 
  mutate(genoAA = case_when(
    copy_num1 >= 1 & copy_num2 >= 1 ~ 1,
    copy_num1 >= 1 & copy_num2 == 0 ~ 0,
    copy_num1 == 0 & copy_num2 >= 1 ~ 0,
    copy_num1 == 0 & copy_num2 == 0 ~ 0
  ),genoAB = case_when(
    copy_num1 >= 1 & copy_num2 >= 1 ~ 0,
    copy_num1 >= 1 & copy_num2 == 0 ~ 1,
    copy_num1 == 0 & copy_num2 >= 1 ~ 1,
    copy_num1 == 0 & copy_num2 == 0 ~ 0
  ),genoBB = case_when(
    copy_num1 >= 1 & copy_num2 >= 1 ~ 0,
    copy_num1 >= 1 & copy_num2 == 0 ~ 0,
    copy_num1 == 0 & copy_num2 >= 1 ~ 0,
    copy_num1 == 0 & copy_num2 == 0 ~ 1)
  )

#summarize PAV genotypes for each gene across all 7 PAV genotype categories
#AAonly, ABonly, BBonly, gAA_AB,gAA_BB,gAB_BB,gAA_AB_BB
genesumtab<-subsetAphHW %>% group_by(Gene,species) %>% reframe(sumAA=sum(genoAA,na.rm = T),sumAB=sum(genoAB,na.rm = T),sumBB=sum(genoBB,na.rm = T)) %>%
  mutate(genoAAonly = case_when(
    sumAA >= 1 & sumAB == 0 & sumBB == 0 ~ 1,
    sumAB == 0 ~ 0,
    sumBB == 0 ~ 0
  ),
  genoABonly = case_when(
    sumAA == 0 & sumAB > 0 & sumBB == 0 ~ 1,
    sumBB > 0 | sumAA > 0 ~ 0
  ),
  genoBBonly = case_when(
    sumBB > 1 & sumAB == 0 & sumAA == 0 ~ 1,
    sumBB == 0 ~ 0
  ),
  gAA_AB = case_when(
    sumAA > 0 & sumAB > 0 & sumBB == 0 ~ 1,
    sumAB == 0 | sumBB == 0 ~ 0
  ),
  gAA_BB = ifelse(
    sumAA > 0 & sumBB > 0 & sumAB == 0, 1, 0
  ),
  gAB_BB = ifelse(
    sumAA == 0 & sumAB > 0 & sumBB > 0, 1, 0
  ),
  gAA_AB_BB = ifelse(
    sumAA > 0 & sumAB > 0 & sumBB > 0, 1,0
  ))

#final summary of CNV genotypes, adding genes across categories as checks

genesumtab %>% filter(species %in% c("AW","AC","AI")) %>%
  group_by(species) %>% reframe(AAonly = sum(genoAAonly,na.rm = T),
                                ABonly = sum(genoABonly,na.rm = T),
                                BBonly = sum(genoBBonly, na.rm = T),
                                ABpres = sum(genoABonly,gAA_AB,gAB_BB,gAA_AB_BB,na.rm = T),
                                BB_pres = sum(genoBBonly,gAA_BB,gAB_BB,gAA_AB_BB,na.rm = T),
                                sum_all_genes=sum(c(genoAAonly,genoABonly,genoBBonly,gAA_AB,gAA_BB,gAB_BB,gAA_AB_BB),na.rm = T))

#now summarize PAVs in 2 or more haplotypes
genesumtab %>% filter(sumAB > 1 | sumBB > 1) %>%
  filter(species %in% c("AW","AC","AI")) %>%
  group_by(species) %>% reframe(AAonly = sum(genoAAonly,na.rm = T),
                                ABonly = sum(genoABonly,na.rm = T),
                                BBonly = sum(genoBBonly, na.rm = T),
                                ABpres = sum(genoABonly,gAA_AB,gAB_BB,gAA_AB_BB,na.rm = T),
                                BB_pres = sum(genoBBonly,gAA_BB,gAB_BB,gAA_AB_BB,na.rm = T),
                                sum_all_genes=sum(c(genoAAonly,genoABonly,genoBBonly,gAA_AB,gAA_BB,gAB_BB,gAA_AB_BB),na.rm = T))


#now summarize PAVs in 6 or more haplotypes
genesumtab %>% filter(sumAB > 5 | sumBB > 5) %>%
  filter(species %in% c("AW","AC","AI")) %>%
  group_by(species) %>% reframe(AAonly = sum(genoAAonly,na.rm = T),
                                ABonly = sum(genoABonly,na.rm = T),
                                BBonly = sum(genoBBonly, na.rm = T),
                                ABpres = sum(genoABonly,gAA_AB,gAB_BB,gAA_AB_BB,na.rm = T),
                                BB_pres = sum(genoBBonly,gAA_BB,gAB_BB,gAA_AB_BB,na.rm = T),
                                sum_all_genes=sum(c(genoAAonly,genoABonly,genoBBonly,gAA_AB,gAA_BB,gAB_BB,gAA_AB_BB),na.rm = T))


#total number of genotypes per species heterozygous for deletion (+1/0)
genesumtab %>% filter(species %in% c("AW","AC","AI")) %>% group_by(species) %>%
  reframe(totBB = sum(sumAB))

#with polyA* filter
genesumtab %>% filter(species %in% c("AW","AC","AI") & sumAB > 1 | sumBB > 1) %>% group_by(species) %>%
  reframe(totBB = sum(sumAB))

#now summarize PAVs in 6 or more haplotypes
genesumtab %>% filter(species %in% c("AW","AC","AI") & sumAB > 5 | sumBB > 5) %>% group_by(species) %>%
  reframe(totBB = sum(sumAB))

#total number of genotypes per species homozygous for deletion (0/0)
genesumtab %>% filter(species %in% c("AW","AC","AI")) %>% group_by(species) %>%
  reframe(totBB = sum(sumBB))

#with polyA* filter
genesumtab %>% filter(species %in% c("AW","AC","AI") & sumAB > 1 | sumBB > 1) %>% group_by(species) %>%
  reframe(totBB = sum(sumBB))

#now summarize PAVs in 6 or more haplotypes
genesumtab %>% filter(species %in% c("AW","AC","AI") & sumAB > 5 | sumBB > 5) %>% group_by(species) %>%
  reframe(totBB = sum(sumBB))
