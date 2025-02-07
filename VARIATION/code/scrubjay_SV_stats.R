library(tidyverse)
library(ggpubr)
library(ggsci)
library("scales")

#Shortened version
svs <- read_table("/n/holylfs05/LABS/informatics/Users/dkhost/scrubjay_analysis/variation_final_v2.sv_short.gte50bp.tab")

#Big version
svs <- read_table("/n/holylfs05/LABS/informatics/Lab/scrubjay/VARIATION/CLEANED_VCFs/variation_final_v2.sv_short.tab")

svs2 <- svs %>%
  mutate(across(c(subtype), ~ifelse(str_detect(subtype, "DEL"), "DEL", .))) %>%
  mutate(across(c(subtype), ~ifelse(str_detect(subtype, "INS"), "INS", .))) %>%
  mutate(sv_length = case_when(base_allele_len >= alt_len_max ~ base_allele_len, base_allele_len < alt_len_max ~ alt_len_max)) %>% #this line creates a new column giving total sv length depending whether ref or alt is longer
  mutate(across(c(subtype), ~ifelse(sv_length < 50, "INDEL (< 50bp)", .)))

  
names(svs2)[names(svs2) == 'repeat'] <- 'repeatID'

#option to remove indels:
svs2 <- svs2 %>% filter(sv_length >= 50)

#test <- svs2 %>% filter(sv_length >= 50)

#Summary table of SV sizes
sizeSum <- svs2 %>%
  group_by(subtype) %>%
  summarise(totalBases = sum(sv_length), totalCounts = n())

write.table(sizeSum, file = '/n/holylfs05/LABS/informatics/Lab/scrubjay/VARIATION/CLEANED_VCFs/PLOTS/svSumTable.tsv', col.names = TRUE, row.names = FALSE, sep = '\t', quote=FALSE)

#SV SIZES HISTOGRAM

p1 <- gghistogram(svs2, x="sv_length", fill="subtype", bins=15, palette="aaas") + 
  facet_wrap(vars(subtype), ncol=3)+
  theme(text = element_text(size=18))+
#  scale_y_continuous(labels = function(x) format(x, scientific = FALSE))+
  scale_y_continuous(breaks=c(0,10,100,1000,10000,100000), labels = function(x) format(x, scientific = FALSE),trans=scales::pseudo_log_trans(base = 10))+
  scale_x_continuous(labels = function(x) format(x/1000, scientific = FALSE))+
  xlab("SV length (kb)")+
  ylab("# SVs")

plot(p1)

ggsave(p1, file="/n/holylfs05/LABS/informatics/Lab/scrubjay/VARIATION/CLEANED_VCFs/PLOTS/scrubjay_SVhisto.pdf", width=11, height = 8.5, device = "pdf")

------

#SV ASSOCIATED WITH GENOMIC REGION
#Bed file of total region start-ends
regionSizes <- read_table("/n/holylfs05/LABS/informatics/Users/dkhost/scrubjay_analysis/genomic_overlaps.merged_d100.bed",col_names = FALSE)
colnames(regionSizes) <- c("chrom","bedStart","bedEnd","region")

regionSizesSum <- regionSizes %>%
  mutate(region = case_when(
    str_detect(region, "cnee") ~ "cnee", #Rather than ignoring SVs with multiple overlaps:
    str_detect(region, "exon") ~ "exon", #In order, see if label contains a given region overlap 
    str_detect(region, "intron") ~ "intron", #And rename that region to it
    str_detect(region, "intergenic") ~ "intergenic",
    str_detect(region, "promoter") ~ "intergenic",
    str_detect(region, "ce") ~ "cnee",
    TRUE ~ "Undefined"
  )) %>%
  group_by(region) %>% summarise(regionTotal = sum(bedEnd-bedStart)) 

regions <- svs2 %>%
  mutate(region = case_when(
    str_detect(overlap, "cnee") ~ "cnee", #Rather than ignoring SVs with multiple overlaps:
    str_detect(overlap, "exon") ~ "exon", #In order, see if label contains a given region overlap 
    str_detect(overlap, "intron") ~ "intron", #And rename that region to it
    str_detect(overlap, "intergenic") ~ "intergenic",
    str_detect(overlap, "promoter") ~ "intergenic",
    str_detect(overlap, "ce") ~ "cnee",
    TRUE ~ "Undefined"
  )) %>%
  group_by(region,subtype) %>% summarise(regionCount = n()) %>% #mutate(count_name_occurr = n())
  left_join(.,regionSizesSum, by = 'region') %>% 
  mutate(svPerMb = (regionCount/regionTotal)*1e6 )

#This plots number of SVs per region type
p2 <- ggplot(regions,aes(x=reorder(region,-regionCount),y=regionCount,fill=region)) + 
  geom_bar(stat='identity',position = 'dodge')+
  facet_wrap(vars(subtype), ncol=3, strip.position = "top")+
  theme_pubr()+
  scale_fill_npg()+
  theme(text = element_text(size=18))+
  theme(axis.text.x=element_text(angle=45,hjust=1,size=18))+
  theme(legend.position="none")+
  #  geom_text(aes(label=count_name_occurr), position=position_dodge(width=0.9), vjust=-0.25)+
  xlab("")+
  ylab("# of SVs")

plot(p2)

ggsave(p2, file="/n/holylfs05/LABS/informatics/Lab/scrubjay/VARIATION/CLEANED_VCFs/PLOTS/scrubjay_SV_perRegion.pdf", width=11, height = 8.5, device = "pdf")

#Output table
write.table(regions, file = '/n/holylfs05/LABS/informatics/Lab/scrubjay/VARIATION/CLEANED_VCFs/PLOTS/SV_perRegionTable.tsv', col.names = TRUE, row.names = FALSE, sep = '\t', quote=FALSE)


#This plots number of SVs per Mb of region type
p3 <- ggplot(regions,aes(x=reorder(region,-svPerMb),y=svPerMb,fill=region)) + 
  geom_bar(stat='identity',position = 'dodge')+
  facet_wrap(vars(subtype), ncol=3, strip.position = "top")+
  theme_pubr()+
  scale_fill_npg()+
  theme(text = element_text(size=18))+
  theme(axis.text.x=element_text(angle=45,hjust=1,size=18))+
  theme(legend.position="none")+
#  geom_text(aes(label=count_name_occurr), position=position_dodge(width=0.9), vjust=-0.25)+
  xlab("")+
  ylab("# of SVs / Mb")

plot(p3)

ggsave(p3, file="/n/holylfs05/LABS/informatics/Lab/scrubjay/VARIATION/CLEANED_VCFs/PLOTS/scrubjay_SV_perRegionMb.pdf", width=11, height = 8.5, device = "pdf")


#This plots the size distribution of SVs within each genomic region
p3v2 <- svs2 %>%
  mutate(region = case_when(
    str_detect(overlap, "cnee") ~ "cnee", #Rather than ignoring SVs with multiple overlaps:
    str_detect(overlap, "exon") ~ "exon", #In order, see if label contains a given region overlap 
    str_detect(overlap, "intron") ~ "intron", #And rename that region to it
    str_detect(overlap, "intergenic") ~ "intergenic",
    str_detect(overlap, "promoter") ~ "intergenic",
    str_detect(overlap, "ce") ~ "cnee",
    TRUE ~ "Undefined"
  )) %>%
  ggplot(.,aes(x=region,y=sv_length,fill=region)) +
  geom_boxplot()+
  facet_wrap(vars(subtype), ncol=3, strip.position = "top")+
  theme_pubr()+
  scale_fill_npg()+
  theme(text = element_text(size=18))+
  theme(axis.text.x=element_text(angle=45,hjust=1,size=18))+
  theme(legend.position="none")+
  scale_y_continuous(breaks=c(0,10,100,1000,10000,100000), labels = function(x) format(x, scientific = FALSE),trans=scales::pseudo_log_trans(base = 10))

plot(p3v2)

ggsave(p3v2, file="/n/holylfs05/LABS/informatics/Lab/scrubjay/VARIATION/CLEANED_VCFs/PLOTS/scrubjay_SV_regionSizeDistrib.pdf", width=11, height = 8.5, device = "pdf")


-------

#REPEAT ASSOCIATED SVs
# repSVs <- svs2 %>% 
#   filter(repeat_family != "none") %>%
#   mutate(across(c(repeat_family), ~ifelse(str_detect(repeat_family, ","), "Multi", .))) %>%
#   mutate(across(c(repeat_family), ~ifelse(str_detect(repeat_family, "DNA"), "DNA", .))) %>%
#   mutate(across(c(repeat_family), ~ifelse(str_detect(repeat_family, "LINE"), "LINE", .))) %>%
#   mutate(across(c(repeat_family), ~ifelse(str_detect(repeat_family, "LTR"), "LTR", .))) %>%
#   mutate(across(c(repeat_family), ~ifelse(str_detect(repeat_family, "SINE"), "SINE", .))) %>%
#   group_by(repeat_family,subtype) %>% mutate(count_name_occurr = n()) %>% #this line groups by family and SV tpye and counts occurrences and adds a column
#   mutate(sv_length = case_when(base_allele_len >= alt_len_max ~ base_allele_len, base_allele_len < alt_len_max ~ alt_len_max)) #this line creates a new column giving total sv length depending whether ref or alt is longer
# 
# test <- svs2[1:24,]

#Read in bed file of repeat annotations to get total length
repeatSizes <- read_table("/n/holylfs05/LABS/informatics/Users/dkhost/scrubjay_analysis/aphWoo.repeats.merged_d10.bed",col_names = FALSE)
colnames(repeatSizes) <- c("chrom","bedStart","bedEnd","repeat_family")

repeatSizesSum <- repeatSizes %>%
  mutate(repeat_family = case_when(
    str_detect(repeat_family, ",") ~ "Multi", #Rather than ignoring SVs with multiple overlaps:
    str_detect(repeat_family, "DNA") ~ "DNA", #In order, see if label contains a given region overlap 
    str_detect(repeat_family, "LTR") ~ "LTR", #And rename that region to it
    str_detect(repeat_family, "SINE") ~ "SINE",
    str_detect(repeat_family, "LINE") ~ "LINE",
    TRUE ~ repeat_family
  )) %>%
  group_by(repeat_family) %>% summarise(repeatTotal = sum(bedEnd-bedStart), repeatTotalCounts = n()) %>%
  filter(repeat_family != 'tRNA' & repeat_family != 'rRNA')


repSVsExpand <- svs2 %>%
  filter(repeat_family != "none") %>%
  separate_longer_delim(c(repeat_family), delim=",") %>%
  distinct() %>% #remove duplicate rows so not double counting SVs with multiple repeats but of same family
  mutate(repeat_family = case_when(
    str_detect(repeat_family, "DNA") ~ "DNA", #In order, see if label contains a given region overlap 
    str_detect(repeat_family, "LTR") ~ "LTR", #And rename that region to it
    str_detect(repeat_family, "SINE") ~ "SINE",
    str_detect(repeat_family, "LINE") ~ "LINE",
    TRUE ~ repeat_family
  )) %>%
  group_by(repeat_family,subtype) %>%
  summarize(numSV = n(), numSVbp = sum(sv_length)) %>%
  left_join(.,repeatSizesSum,by='repeat_family') %>%
  filter(!grepl("RNA",repeat_family)) %>% 
  mutate(SVperFamMb = (numSV/repeatTotal)*1e6, SVbyFamCounts = numSV/repeatTotalCounts)


# repSVsSum <- repSVs %>% group_by(repeat_family,subtype) %>%
#   filter(repeat_family != 'tRNA' & repeat_family != 'rRNA') %>%
#   summarise(sv_length = sum(sv_length), familyCount = n()) %>%
#   left_join(.,repeatSizesSum, by='repeat_family') %>%
#   mutate(svPerMb = (familyCount/repeatTotal)*1e6 ) %>%
#   mutate(svLengthPerMb = (sv_length/repeatTotal)*1e6 )
  

#This plots the number of bases per repeat-associated SV type
p4 <- ggplot(repSVsExpand,aes(x=reorder(repeat_family,-numSV),y=numSV,fill=repeat_family))+
  geom_bar(stat='identity')+
  theme_pubr()+
  scale_fill_lancet()+
#  theme(axis.text.x=element_text(angle=45,hjust=1,size=14))+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  theme(legend.position="right")+
  theme(text = element_text(size=18))+
  facet_grid(cols = vars(subtype))+
#  scale_y_continuous(labels = function(x) format(x/1e6, scientific = FALSE))+
#  geom_text(aes(label=paste(format((sv_length/1e6), digits = 2), "Mb")), position=position_dodge(width=0.9), vjust=-0.25)+
  xlab("")+
  ylab("# of SVs")

plot(p4)

ggsave(p4, file="/n/holylfs05/LABS/informatics/Lab/scrubjay/VARIATION/CLEANED_VCFs/PLOTS/scrubjay_SV_perRepeatFam.pdf", width=11, height = 8.5, device = "pdf")

#Output table
write.table(repSVsExpand, file = '/n/holylfs05/LABS/informatics/Lab/scrubjay/VARIATION/CLEANED_VCFs/PLOTS/SV_perRepeatFamTable.tsv', col.names = TRUE, row.names = FALSE, sep = '\t', quote=FALSE)


#This plots number of SVs per MB of repeat family type
p5 <- ggplot(repSVsExpand,aes(x=reorder(repeat_family,-SVperFamMb),y=SVperFamMb,fill=repeat_family))+
  geom_bar(stat='identity')+
  theme_pubr()+
  scale_fill_lancet()+
#  theme(axis.text.x=element_text(angle=45,hjust=1,size=24))+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  theme(legend.position="right")+
  theme(text = element_text(size=18))+
  facet_grid(cols = vars(subtype))+
 # scale_y_continuous(labels = function(x) format(x, scientific = FALSE))+
  #  geom_text(aes(label=paste(format((sv_length/1e6), digits = 2), "Mb")), position=position_dodge(width=0.9), vjust=-0.25)+
  xlab("")+
  ylab("# SVs / Repeat family Mb")

plot(p5)

ggsave(p5, file="/n/holylfs05/LABS/informatics/Lab/scrubjay/VARIATION/CLEANED_VCFs/PLOTS/scrubjay_SV_perRepeatFamMb.pdf", width=11, height = 8.5, device = "pdf")


#This plots # SV divided by # of repeat family annotations
p5v2 <- ggplot(repSVsExpand,aes(x=reorder(repeat_family,-SVbyFamCounts),y=SVbyFamCounts,fill=repeat_family))+
  geom_bar(stat='identity')+
  theme_pubr()+
  scale_fill_lancet()+
  #  theme(axis.text.x=element_text(angle=45,hjust=1,size=24))+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  theme(legend.position="right")+
  theme(text = element_text(size=18))+
  facet_grid(cols = vars(subtype))+
  # scale_y_continuous(labels = function(x) format(x, scientific = FALSE))+
  #  geom_text(aes(label=paste(format((sv_length/1e6), digits = 2), "Mb")), position=position_dodge(width=0.9), vjust=-0.25)+
  xlab("")+
  ylab("# SVs / # repeat family")

plot(p5v2)

ggsave(p5v2, file="/n/holylfs05/LABS/informatics/Lab/scrubjay/VARIATION/CLEANED_VCFs/PLOTS/scrubjay_SV_byRepeatFamCount.pdf", width=11, height = 8.5, device = "pdf")


#Looking at repeats per region, split by repeat type
SVbyRegionRepeats <- svs2 %>%
  separate_longer_delim(c(repeat_family), delim=",") %>%
  distinct() %>%
  mutate(region = case_when(
    str_detect(overlap, "cnee") ~ "cnee", #Rather than ignoring SVs with multiple overlaps:
    str_detect(overlap, "exon") ~ "exon", #In order, see if label contains a given region overlap 
    str_detect(overlap, "intron") ~ "intron", #And rename that region to it
    str_detect(overlap, "intergenic") ~ "intergenic",
    str_detect(overlap, "promoter") ~ "intergenic",
    str_detect(overlap, "ce") ~ "cnee",
    TRUE ~ "Undefined"
  )) %>%
  mutate(repeat_family = case_when(
    str_detect(repeat_family, "DNA") ~ "DNA", #In order, see if label contains a given region overlap 
    str_detect(repeat_family, "LTR") ~ "LTR", #And rename that region to it
    str_detect(repeat_family, "SINE") ~ "SINE",
    str_detect(repeat_family, "LINE") ~ "LINE",
    TRUE ~ repeat_family
  )) %>%
  filter(!grepl("RNA",repeat_family)) %>%
  group_by(region,repeat_family,subtype) %>% summarise(regionCount = n()) %>% #mutate(count_name_occurr = n())
  left_join(.,regionSizesSum, by = 'region') %>%
  mutate(SVperRegionMb = (regionCount/regionTotal)*1e6)


show_col(pal_lancet()(9))
pal_lancet()(9)
lancet_custom <- c("#00468BFF","#ED0000FF","#42B540FF","#0099B4FF","white","#925E9FFF","#FDAF91FF","#AD002AFF","#ADB6B6FF","#1B1919FF")


p6 <- ggplot(SVbyRegionRepeats,aes(x=reorder(region,-SVperRegionMb),y=SVperRegionMb,fill=repeat_family))+
  geom_bar(stat='identity',position='stack',color="black")+
  theme_pubr()+
  #scale_fill_lancet()+
  scale_fill_manual(values = lancet_custom)+
  # theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  theme(legend.position="right")+
  theme(text = element_text(size=18))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  facet_grid(cols = vars(subtype))+
  # scale_y_continuous(labels = function(x) format(x, scientific = FALSE))+
  #  geom_text(aes(label=paste(format((sv_length/1e6), digits = 2), "Mb")), position=position_dodge(width=0.9), vjust=-0.25)+
  xlab("")+
  ylab("# SVs / Mb")

plot(p6)

#Output table
write.table(SVbyRegionRepeats, file = '/n/holylfs05/LABS/informatics/Lab/scrubjay/VARIATION/CLEANED_VCFs/PLOTS/SV_perRegion_byRepeatTable.tsv', col.names = TRUE, row.names = FALSE, sep = '\t', quote=FALSE)

ggsave(p6, file="/n/holylfs05/LABS/informatics/Lab/scrubjay/VARIATION/CLEANED_VCFs/PLOTS/scrubjay_SV_perRegion_byRepeat.pdf", width=11, height = 8.5, device = "pdf")
