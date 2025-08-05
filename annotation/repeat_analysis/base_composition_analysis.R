library(ggpubr)
#example read in of output from seqtk comp to sum by individual or haplotype and format for plotting
#the file sat_sj_circ1_2268_comp_files_all_haps.txt is a concatenation of seqtk comp outputs for each haplotype specifically for the common satellite sj_circ1_2268
#the file all_haps_seqtk_basecomp.txt is a concatenation of seqtk comp outputs for each haplotype, genome wide

bc<-read.table("sat_sj_circ1_2268_comp_files_all_haps.txt",sep = "\t",comment.char = "@")
allcomp<-read.table("all_haps_seqtk_basecomp.txt",sep = "\t",comment.char = "@")

#format both files to eventually merge
bc$species<-str_before_first(bc$V1,"_")
allcomp$species<-str_before_first(allcomp$V1,"_")
allcomp$comp<-"all_haps"
bc$comp<-"circ1_2268"
bc$len<-bc$V3-bc$V2
bc<-bc[,c(1,17,4:7,16)]

allcomp<-allcomp[,c(1:6,15)]
colnames(allcomp)<-c("contig","len","#A","#C","#G","#T","comp")
colnames(bc)<-c("contig","len","#A","#C","#G","#T","comp")

#merge files and format further
comp<-rbind(bc,allcomp)
comp$species<-str_before_first(comp$contig,"_")
comp$hap<-gsub(".p_ctg","",str_before_nth(comp$contig,"#",2))
comp$indiv<-str_before_first(comp$hap,"#")

#comp is an example of a file summarizing base composition across the whole genome for each individual or haplotype


#construction of base composition table for plots and statistical analysis

#this file contains base composition counts for all satellite regions combined
allsatbc<-read.table("sum_AGCT_seqtk_comp_all_sats_by_indiv.txt",sep = "\t",header = T)

#create base composition files of each of the top five most common satellites individually
#the input files here are simply derived from the comprehensive base composition files available on Dryad
a<-read.table("seqtk_comp_all_sum_haps_sj_sat_circ1_2268.txt",sep = "\t",comment.char = "@",header = T)
a$component<-"sj_sat_circ1_2268"
b<-read.table("seqtk_comp_all_sum_haps_sj_sat_circ2_2119.txt",sep = "\t",header = T)
b$component<-"sj_sat_circ2_2119"
c<-read.table("seqtk_comp_all_sum_haps_sj_sat_circ3-2178.txt",sep = "\t",header = T)
c$component<-"sj_sat_circ3-2178"
d<-read.table("seqtk_comp_all_sum_haps_sj_sat_circ30−18193.txt",sep = "\t",header = T)
d$component<-"sj_sat_circ30−18193"
e<-read.table("seqtk_comp_all_sum_haps_sj_sat_circ4_2216.txt",sep = "\t",header = T)
e$component<-"sj_sat_circ4_2216"
f<-read.table("seqtk_comp_all_sum_haps_sj_sat_circ7-2218.txt",sep = "\t",header = T)
f$component<-"sj_sat_circ7-2218"
g<-read.table("seqtk_comp_sums_all_haps_whole_genome_haps.txt",sep = "\t",header = T)
g$component<-"genome"

allsatbc$component<-"all_sats"
allgc<-rbind(a,b,c,d,e,f,g,allsatbc)

write.table(allgc,"five_sats_whole_genome_seqtk_comp_sums.txt",sep = "\t",col.names = T,row.names = F,quote = F)

#this file contains base composition counts for all non-satellite regions of each haplotype
nonsat<-read.table("seqtk_comp_all_sum_haps_non_sats.txt",sep = "\t",comment.char = "@",header = T)
allgc<-read.table("five_sats_whole_genome_seqtk_comp_sums.txt",sep = "\t",header = T)
allgc <- allgc %>% mutate(GC=(sumC+sumG)/sumlen)
nonsat<-nonsat %>% mutate(GC=(sumC+sumG)/sumlen)
allgc$species<-str_before_first(allgc$indiv,"_")
nonsat$species<-str_before_first(nonsat$indiv,"_")
nonsat$component<-"nonsat"
nonsat<-nonsat[,c(1:6,9,7,8)]

allgc<-rbind(allgc,nonsat)
write.table(allgc,"five_sats_whole_genome_seqtk_comp_sums.txt",sep = "\t",col.names = T,row.names = F,quote = F)

#some exploratory plots
#six panels of GC content for each species (AA, AC, AI, AW, CS and CY)
allgc %>% ggplot(aes(x=component,y=GC,color=component)) + geom_jitter(size=0.8,alpha=0.5) +
  facet_wrap(.~species,ncol=3,scales = "free_y") + theme_bw() +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.95,angle = 90,size = 10))
ggsave("five_sats_genome_GC_content_seqtk_comp_facet_by_species.pdf",height = 6,width = 8)


#plot for AW individually
allgc %>% filter(species=="AW") %>% ggplot(aes(x=component,y=GC,color=component)) + 
  geom_jitter() +theme(axis.text.x = element_text(vjust = 0.5, hjust=0.95,angle = 90,size = 10))

#examine GC content for each genomic component for each species
allgc %>% ggplot(aes(x=species,y=GC,color=species)) + 
  geom_jitter(size=0.8,alpha=0.5) +theme(axis.text.x = element_text(vjust = 0.5, hjust=0.95,angle = 90,size = 10)) +
  facet_wrap(.~component,ncol=3,scales = "free_y") + theme_bw()
ggsave("five_sats_genome_GC_content_seqtk_comp_facet.pdf",height = 6,width = 8)

#prepare wider GC content table for plotting
allgct<-allgc %>% select(c(1,6:9))
allgc<-allgc %>% mutate(GCcount=sumC+sumG)
allgct<-allgc %>% select(c(1,6:10))
allgctt<-allgct %>% pivot_wider(names_from = component,values_from = c(sumlen,GCcount,GC))

write.table(allgctt,"five_sats_non_sat_genome_GC_content_seqtk_comp_wider.txt",sep = "\t",col.names = T,quote = F,row.names = F)
allgctt<-read.table("five_sats_non_sat_genome_GC_content_seqtk_comp_wider.txt",header = T,sep = "\t")


`%nin%` = Negate(`%in%`)
colorscheme<-c("#114EF2",brewer.pal(5,"RdYlBu")[1],"#8AA5EE","#A9ADB8","lightgreen","black")
p<-allgctt %>% ggplot(aes(x=sumlen_all_sats/1e6,y=GC_genome,color=species)) + geom_point() +
  ylab("GC content per individual") + xlab("total Mbp satellite per individual") +
  theme_bw() + scale_color_manual(values = colorscheme)
ggsave("scatterplot_all_hap_basecomp_vs_tot_sat_abun_2.pdf",width = 6,height = 6)

allgctt$species<-factor(allgctt$species,levels = c("AI","AA","AW","AC","CS","CY"))

q<-allgctt %>% ggplot(aes(y=GC_genome,x=species,fill=species)) + geom_boxplot() +
  ylab("GC content per individual") + xlab("species") +
  scale_fill_manual(values = colorscheme) +
  theme_bw() + theme(legend.position="none")
ggsave("boxplot_all_hap_basecomp.pdf",width = 6,height = 6)

values = rev(brewer.pal(11,name="RdYlBu")[c(1,2,11,3,7)])
r<-allgctt %>% ggplot(aes(y=GC_all_sats,x=species,fill=species)) + geom_boxplot() +
  ylab("GC content of satellite component per individual") + xlab("species") +
  scale_fill_manual(values = colorscheme) +
  theme_bw() + theme(legend.position="none")
ggsave("boxplot_all_hap_sats_only_basecomp.pdf",width = 6,height = 6)

s<-allgc %>% filter(component %nin% c("all_sats","genome")) %>% ggplot(aes(x=species,y=GC,color=species)) + 
  geom_jitter(size=0.8,alpha=0.5) +theme(axis.text.x = element_text(vjust = 0.5, hjust=0.95,angle = 90,size = 10)) +
  facet_wrap(.~component,ncol=2,scales = "free_y") + theme_bw()# + theme(legend.position="none") 
ggarrange(q,s,r,p,
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2)
ggsave("four_panels_sideways_basecomp_fin_alt3.pdf",width = 12,height = 10)

