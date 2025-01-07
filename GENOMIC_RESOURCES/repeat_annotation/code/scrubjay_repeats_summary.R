library(tidyverse)
library(ggpubr)
library(ggsci)
library(pheatmap)
library(RColorBrewer)

ref <- read_table("/n/holylfs05/LABS/informatics/Lab/scrubjay/GENOMIC_RESOURCES/repeat_annotation/results/aphWoo.repeats.tab",col_names = FALSE)

refFamSum <- ref %>%
  rename("chrom" = "X1", "bedStart" = X2, "bedEnd" = X3, "repID" = "X4", "repFamily" = "X5") %>%
  mutate(repFamily = case_when(
    str_detect(repFamily, "DNA") ~ "DNA", #In order, see if label contains a given region overlap 
    str_detect(repFamily, "LTR") ~ "LTR", #And rename that region to it
    str_detect(repFamily, "SINE") ~ "SINE",
    str_detect(repFamily, "LINE") ~ "LINE",
    TRUE ~ repFamily
  )) %>%
  group_by(repFamily) %>%
  summarise(famRepCount = n(), famRepBP = sum(bedEnd-bedStart))

lancet_custom <- c("#00468BFF","#ED0000FF","#42B540FF","#0099B4FF","#925E9FFF","#FF66CC","#FDAF91FF","#AD002AFF","#ADB6B6FF","#FFFFCC","#1B1919FF")

ggplot(refFamSum, aes(x=reorder(repFamily,-famRepCount),y=famRepCount,fill=repFamily)) +
  geom_bar(stat='identity')+
  theme_pubr()+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_fill_manual(values = lancet_custom)

p1 <- ggplot(refFamSum, aes(x=reorder(repFamily,-famRepBP),y=famRepBP,fill=repFamily)) +
  geom_bar(stat='identity',color='black')+
  theme_classic2()+
  theme(text = element_text(size=18))+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
 # scale_y_continuous(breaks=c(0,1e7,1e8), labels = function(x) format(x, scientific = FALSE))+
  scale_y_continuous(breaks=c(0,25e6,50e6,75e6,100e6,125e6), labels = function(x) format(x/1e6, scientific = FALSE))+
  ylab("Mb per family")+
  scale_fill_manual(values = lancet_custom)
plot(p1)

ggsave(p1,file="/n/holylfs05/LABS/informatics/Lab/scrubjay/GENOMIC_RESOURCES/repeat_annotation/results/scrubjay_mbPerRepeatFamily_bar.pdf", width=11, height = 8.5, device = "pdf")
