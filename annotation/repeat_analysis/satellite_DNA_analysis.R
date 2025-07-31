#access *.len files from Satellite Repeat Finder analysis of fastq files and extract relevant columns

setwd("/Users/scott/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Scrub\ Jays/satellites/Heng_analysis/March_2024_analyses/sj_sats-abun-20240823")
library(strex)
library(dplyr)
library(tidyverse)
i<-1
satreads<-data.frame()
for (i in 1:length(dir(pattern=".len"))){
  x<-read.table(dir(pattern=".len")[i],comment.char = "@")
  x<-x[,c(1,2,4,5)]
  x$indiv<-NA
  x$indiv<-str_before_last_dot(dir()[i])
  x$unit_length<-NA
  x$unit_length<-as.numeric(str_after_last(x$V1,"-"))
  satreads<-rbind(satreads,x)
}

satreads$species<-str_before_first(satreads$indiv,"_")
satreads<-satreads %>% mutate(length_category = case_when(unit_length < 100 ~ "< 100",
                                                      unit_length > 100 & unit_length < 1000 ~ "100-1000",
                                                      unit_length > 1000 & unit_length < 2000 ~ "1000-2000",
                                                      unit_length > 2000 & unit_length < 5000 ~ "2000-5000",
                                                      unit_length > 5000 & unit_length < 10000 ~ "5000-10000",
                                                      unit_length > 10000 ~ ">10000"))


colnames(satreads)[1:4]<-c("satellite","tot_bp_reads","without_long","with_long")
write.table(satreads,"Heng_summary_satellite_table_Aug_2024.txt",sep = "\t",quote = F,col.names = T,row.names = F)

#combine info from SRF tables with primary and haplotype assembly lengths and with SRF results from assemblies
satreads<-read.table("/Users/scott/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Scrub Jays/satellites/Heng_analysis/March_2024_analyses/sj_sats-abun-20240823/Heng_summary_satellite_table_Aug_2024.txt",header = T,sep = "\t",comment.char = "@")
satass<-read.table("/Users/scott/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Scrub Jays/satellites/Danielle_srf_files_2/final_contig_srf_Danielle/sj_sats_genomes_CONCAT.len",sep = "\t",comment.char = "@")
haps<-read.table("/Users/scott/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Scrub\ Jays/assemblies/danielle_assembly_stats/ASSEMBLY_STATS.tsv",header = T)

colnames(satass)[c(1,2,6)]<-c("sat","tot_bp_in_assem","indiv")
satass<-satass[,c(1,2,6)]

#read in list of 28 most common satellites
commsats<-read.table("/Users/scott/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Scrub Jays/satellites/top_28_sat_list_from_heatmap.txt",comment.char = "@")
commsats$V1<-gsub("sat_","sat#",commsats$V1)
satass <- satass %>% filter(sat %in% commsats$V1)
satass$species<-str_before_first(satass$indiv,"_")

#check concordance of sample names
setdiff(satass$indiv,haps$Sample)
setdiff(haps$Sample,satass$indiv)
length(unique(satass$sat))

#fold in primary assembly lengths for each individual
for(i in 1:length(satass$indiv)){
  matchnum<-which(haps$Sample==satass$indiv[i])
  satass$pri_assem_len[i]<-haps$pri.Tot_len[matchnum]
}

colnames(satass)[c(1:3)]<-c("satellite","tot_bp_assem","indiv")

satreads<-satreads %>% filter(satellite %in% commsats$V1)

#add satellite proportion in assemblies ('tot_bp_assem') to satreads
satreads$tot_bp_assem<-NA
for(i in 1:dim(satreads)[1]){
  ifelse(dim(match_df(satass[,c(1,3)],satreads[i,c(1,5)]))[1] == 0,i<-i+1,
         pick<-as.numeric(rownames(match_df(satass[,c(1,3)],satreads[i,c(1,5)]))))
  satreads$tot_bp_assem[i]<-satass$tot_bp_assem[pick]
}

#rearrange and streamline satread table
colnames(satreads)
satreads2<-satreads[,c(1,6,5,7,2,4,9,10,8)]
write.table(satreads2,"final_satellite_summary_reads_assemblies_dryad.txt",sep = "\t",col.names = T,row.names = F,quote = F)

#make supplementary figure S13
satreads2 %>% ggplot(aes(x=with_long,y=tot_bp_assem/pri_assem_len,label=satellite,fill=species,color=species)) + geom_point(aes(shape = satellite), size = 2) + 
  ylab("proportion of satellites in primary assembly") +
  xlab("proportion of satellites in raw reads")  +
  scale_shape_manual(values=c(15,rep(19,3),18,19,19,17,rep(19,8),3,4,rep(19,10)),
                     guide = guide_legend(override.aes = list(size = 3))) +
  scale_color_manual(values = c("#A9ADB8","#114EF2","#8AA5EE","black")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),panel.background = element_blank()) +
  geom_abline(slope = 1,color="gray",linetype="dashed") + theme(legend.position="none")
ggsave("/Users/scott/Library/CloudStorage/OneDrive-HarvardUniversity/Research/Scrub\ Jays/satellites/final_reads_vs_assems_top_28_sats_shapes_line_noleg_v2.pdf",p,height = 8,width = 8)
