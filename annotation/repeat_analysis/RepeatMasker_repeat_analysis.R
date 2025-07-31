library(dplyr)
library(strex)
library(tidyverse)
library(tidyr)
library(RColorBrewer)

#code used for making supplementary table S4 and associated repeat plots
#start with concatenation of RepeatMasker output bedfiles, placed on Dryad

#script used to prepare final repmask files for Zenodo July 2025
setwd("/n/netscratch/edwards_lab/Lab/sedwards/scrubjay/repmask/repmask_hap_nornd_bedfiles")
library(strex)
library(dplyr)
library(tidyverse)
length(dir(pattern = "out.bed"))
for(i in 1:length(dir(pattern = "out.bed"))){
  x<-read.table(dir(pattern = "out.bed")[i],sep = "\t",comment.char = "@",fill = TRUE,row.names=NULL)
  x<-x[,c(1:4,6)]
  colnames(x)<-c("contig","start","end","rep_type","repeat_type")
  x$hap<-str_before_nth(x$contig,"#",2)
  x$species<-str_before_first(x$contig,"_")
  write.table(x,"all_haps_repmask_nornd_with_CS_CY_SJ.bed",sep = "\t",col.names = F,row.names = F,quote = F,append = T)
  print(i)
}

#read in file again and pare down for storage and analysis
reps<-read.table("all_haps_repmask_nornd_with_CS_CY_SJ.bed",sep = "\t",comment.char = "@")
colnames(reps)<-c("contig","start","end","rep_type","repeat_type","hap","species")
reps2<-reps %>% mutate(replen=end-start) %>% group_by(hap,repeat_type) %>% reframe(hap=hap,repeat_type=repeat_type,totreplen=sum(replen),species=species) %>% distinct()
write.table(reps2,"all_haps_repmask_nornd_sums.bed",sep = "\t",col.names = T,row.names = F,quote = F)

#the file all_haps_repmask_nornd_sums.bed was then gzipped and uploaded to Dryad as all_haps_repmask_nornd_sums.bed.gz

#summarize concatenated RepeatMasker bedfile into individual haplotypes
  
samp2<-read.table("/Users/scott/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Scrub\ Jays/RepMask/repmask_all_sj/repmask_all_haps_nornd/all_haps_repmask_nornd_sums_July_2025.bed",sep = "\t",comment.char = "@",header = T)
#need to input major repeat types

#tidy up a bit and get order of listing of all repeat sub-types
unique(samp2$hap)
head(samp2)
tail(samp2)
samp2$hap[which(samp2$hap == "AA_SRR23446543#hap1.p_ctg")]<-"AA_SRR23446543#hap1"
samp2$hap[which(samp2$hap == "AA_SRR23446543#hap2.p_ctg")]<-"AA_SRR23446543#hap2"
samp2$maj_rep_type<-NA

#input major repeat types
samp2[which(samp2$repeat_type %in% c(unique(samp2$repeat_type)[1:8])),5]<-"DNA"
samp2[which(samp2$repeat_type %in% c(unique(samp2$repeat_type)[9:11])),5]<-"LINEs"
samp2[which(samp2$repeat_type %in% c(unique(samp2$repeat_type)[12:19])),5]<-"LTR"
samp2[which(samp2$repeat_type %in% c(unique(samp2$repeat_type)[20])),5]<-"Low_complexity"
samp2[which(samp2$repeat_type %in% c(unique(samp2$repeat_type)[21])),5]<-"RC"
samp2[which(samp2$repeat_type %in% c(unique(samp2$repeat_type)[22:25])),5]<-"SINEs"
samp2[which(samp2$repeat_type %in% c(unique(samp2$repeat_type)[26])),5]<-"Satellite"
samp2[which(samp2$repeat_type %in% c(unique(samp2$repeat_type)[27])),5]<-"Simple"
samp2[which(samp2$repeat_type %in% c(unique(samp2$repeat_type)[28])),5]<-"Unknown"
samp2[which(samp2$repeat_type %in% c(unique(samp2$repeat_type)[29:30])),5]<-"smallrna"

#make barplot from table summary
reptab<-samp2 %>% group_by(maj_rep_type,hap) %>% dplyr::reframe(sumrep=sum(totreplen),hap=hap,species=species) %>% distinct()
reptab %>% ggplot(aes(x=hap,y=sumrep/1e6,fill=maj_rep_type)) + geom_bar(stat="identity",position="stack") +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.95,angle = 90,size = 5)) +
  ylab("megabases") + scale_y_continuous(limits = c(0,400), expand = c(0, 0)) +
  scale_fill_manual(values = rev(brewer.pal(11,name="RdYlBu")[c(1,11,2,10,3,9,4,8,5,7)])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),panel.background = element_blank())

#input individuals for each haplotype
samp2$indiv<-str_before_first(samp2$hap,"#")

#input lengths of each haplotype assembly
haps<-read.table("/Users/scott/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Scrub\ Jays/assemblies/danielle_assembly_stats/ASSEMBLY_STATS.tsv",header = T)
i<-1
haps2<-haps %>% select(c(1,7,10)) %>% pivot_longer(c(2,3),names_to = "haplotype",values_to = "haplen")
haps2$hapname<-gsub(".Tot_len","",paste(haps2$Sample,haps2$haplotype,sep="#"))

#check haplotype names match
setdiff(haps2$hapname,samp2$hap)
setdiff(samp2$hap,haps2$hapname)

#create column for haplotype assembly lengths
samp2$hap_assem_len<-NA

#add haplotype assembly lengths to appropriate rows
for(i in 1:dim(haps2)[1]){
  matchnum<-which(samp2$hap==haps2$hapname[i])
  samp2$hap_assem_len[matchnum]<-haps2$haplen[i]
}

write.table(samp2,"/Users/scott/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Scrub Jays/RepMask/repmask_all_sj/repmask_all_haps_nornd/all_haps_samp2_repmask_nornd_sums.txt",sep = "\t",col.names = T,row.names = F,quote = F)

samp3 <- samp2 %>% group_by(hap,maj_rep_type) %>% reframe(sum_bp = sum(totreplen),hap_assem_len = hap_assem_len,species=species) %>% distinct()
haplevs<-c(31:60,1:2,61:90,3:30,91:94)
#stacked barplot, by haplotype
#Fig2A, main text
pdf("barplot_stacked_major_repeat_abundances_by_haplotype_arranged.pdf", width = 8, height = 6)
samp3 %>% filter(species != "SJ") %>% ggplot(aes(x=factor(hap,levels = unique(hap)[haplevs]),y=sum_bp/1e6,fill=maj_rep_type)) + geom_bar(stat="identity",position="stack") +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.95,angle = 90,size = 5)) +
  ylab("megabases") + scale_y_continuous(limits = c(0,400), expand = c(0, 0)) +
  scale_fill_manual(values = rev(brewer.pal(11,name="RdYlBu")[c(1,11,2,10,3,9,4,8,5,7)])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),panel.background = element_blank())
dev.off()

#facetted barplot, by haplotype
#Supplementary Fig. S11
legend_title<-"major repeat types"
pdf("barplot_facetted_major_repeat_abundances_by_haplotype_arranged.pdf", width = 6, height = 8)
samp3 %>% filter(species != "SJ") %>% ggplot(aes(x=factor(hap,levels = unique(hap)[haplevs]),y=sum_bp/1e6,fill=maj_rep_type)) + geom_bar(stat="identity") +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.95,angle = 90,size = 1.45)) +
  ylab("megabases per haplotype") + xlab("haplotype") + 
  scale_fill_manual(legend_title,values = rev(brewer.pal(11,name="RdYlBu")[c(1,11,2,10,3,9,4,8,5,7)])) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),panel.background = element_blank()) +
  facet_wrap(.~maj_rep_type,ncol=2,scales = "free_y") + theme(legend.position="none")
dev.off()

#bar plot of length of single copy regions by haplotype
#Supplementary Fig. S10
samp3 %>% filter(species != "SJ") %>% group_by(hap) %>% reframe(species=species,single_copy=hap_assem_len-sum(sum_bp)) %>% distinct() %>%
  ggplot(aes(x=hap,y=single_copy,fill=species)) + geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.95,angle = 90,size = 3))

