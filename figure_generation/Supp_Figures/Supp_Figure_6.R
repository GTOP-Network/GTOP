#==============================================#
# Small variant /SV /TR #
# Supp-Figure-6#
#==============================================#
library(data.table)
library(ggplot2)
library(stringi)
library(stringr)
library(dplyr)
library(ggsci)
library(tidyverse)
library(ggpubr)
library(magrittr)
library(scales)
setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig2567/input")
source("geom_boxplot2.R")

# Supp.Fig.6a-c:  Small variant /SV /TR number in diff LRS coverage ------------------------------------------------------------------

#Supp.Fig.6a:  Small variant
number <- fread("Supp_Fig6a.small_variant_downsample.txt")
number$coverage <- factor(number$coverage,levels = c("5x","10x","15x","20x","25x","30x","35x","40x"))
result <- number %>% pivot_longer(cols = c(overlap_num, novel_num),names_to = "category",values_to = "count") 
result$type <- factor(result$type,levels = c("SNP","INDEL"))
result <- result %>% group_by(coverage, type) %>%
  mutate(percentage = count / total_num * 100)
p1 <- ggplot(result, aes(x = type, y = count, fill = category)) + 
  geom_bar(stat="identity", color = "black", width = 1) + 
  theme_classic() + 
  theme(axis.text = element_text(color = "black", size = 12), 
        axis.ticks = element_line(color = "black"),
        legend.position = "top") + 
  labs(x = "LRS Coverage", y = "LRS SNV/INDEL number")+
  facet_wrap(~coverage, ncol = length(unique(result$coverage))) 

#Extended Fig.6b:  TR
number <- fread("Supp_Fig6c.TR_downsample.txt",header = F)
number$V3 <- factor(number$V3,levels = c("5x","10x","15x","20x","25x","30x","35x","40x"))
ggplot(number, aes(x = V3, y = V2)) +
  geom_point() +
  geom_line(group = 1) + 
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(color = "black", size = 12), 
    axis.ticks = element_line(color = "black")) +
  labs(x = "LRS Coverage", y = "TR number")

#Extended Fig.6c:  SV，
ggplot(result[result$type!="SNP",], aes(x = type, y = count, fill = category)) + 
  geom_bar(stat="identity", color = "black", width = 1) + 
  theme_classic() + 
  theme(axis.text = element_text(color = "black", size = 12), 
        axis.ticks = element_line(color = "black"),
        legend.position = "top") + 
  facet_wrap(~coverage, ncol = length(unique(result$coverage))) 

# Extended Fig.6d: small variant genotype accuarcy between LRS and SRS  --------------------------------------------------------

# Extended Fig.6e: SV compare in sniffles2/pbsv/cuteSV  --------------------------------------------------------


# Extended Fig.6f: SV compare with reported datasets  --------------------------------------------------------

dat <- fread("Supp_Fig6f.SV_comapre_with_Other_Datasets.txt")
dat$overlap_category <- str_replace(dat$overlap_category,"overlap_","")
dat <- dat %>% mutate(overlap_category=factor(overlap_category,levels=c("GTEX","gnomAD",
                                                           "1KG","GMTiP_SRS",
                                                           "Han_Chinese")),
         variant_type=factor(variant_type,levels=rev(c("Reported","INS","DEL","BND","DUP","INV"))))
ggplot(dat, aes(x = proportion, y = overlap_category, fill = variant_type)) +
  geom_col(position = "stack") +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Proportion") +
  theme_classic() +
  scale_fill_manual(values = c("DEL" = '#5C87A6', "INS" = '#C56364', "DUP" = '#EEA3A3',
                               "INV" = "#BDDEF2", "BND" = "#80BDE9", "Reported" = "grey"))+
  theme(axis.line=element_line(color='black'),legend.position = "top",
        axis.text=element_text(color='black',size=12))


# Extended Fig.6g-i: LRS SV mutation pattern ---------------------------------------------------

dat <- fread("Supp_Fig6e-g.LRS_SV_info.txt",header = T)
dat <- dat %>% setnames(c("chr","start","ID","length","maf","hwe","end")) %>%
  mutate(svtype=word(ID,start=3,end=3,sep="\\_"))

#Extended Fig.6g: maf
ggplot(dat, aes(maf)) +
  geom_histogram(binwidth = 0.01, fill = "#A2B5CD", color = "black",alpha=0.8) + 
  labs(x = "Minor allele frequency (log scale)", y = "# of SVs") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 12), 
        axis.ticks = element_line(color = "black"))

#Extended Fig.6h: length
dat$length <- abs(as.numeric(dat$length))
p_binwidth <- as.numeric(0.05)
dat$length <- log(abs(dat$length),10)
types <- c('DEL','INS','DUP','INV')
dat$type <- factor(dat$svtype,levels=types,ordered=TRUE)
colors_0 <- c("DEL" = '#5C87A6', "INS" = '#C56364', "DUP" = '#EEA3A3',"INV" = "#BDDEF2", "BND" = "#80BDE9")
ggplot(dat,aes(x=length,group=type)) +
  geom_histogram(aes(fill=type, color=type),binwidth=p_binwidth) +
  xlab('SV Length (bp)') + ylab('# of SVs') +
  scale_fill_manual(breaks=types,values=colors_0) +
  scale_color_manual(breaks=types,values=colors_0) +
  scale_x_continuous(breaks=c(2,3,4,5),labels=c(expression('10'^2),expression('10'^3),expression('10'^4),expression('10'^5))) +
  geom_vline(xintercept= log(c(300,2500,6000),10),linetype='dashed') +
  coord_cartesian(xlim=c(1,6))+
  theme_classic() +
  theme(axis.line=element_line(color='black'),
        axis.text=element_text(color='black',size=rel(1.25)))

#Extended Fig.6i: pergenome count


datt<-fread("Supp_Fig5k.txt")

ggplot(datt, aes(x = sampleID,  y = sv_type_number,fill=type)) +
  geom_col(position = "stack", width = 1) +  
  labs(x = "Sample", y = "# of SVs") +
  theme_classic() +
  scale_fill_manual(values = c("DEL" = '#5C87A6', "INS" = '#C56364', "DUP" = '#EEA3A3',
                               "INV" = "#BDDEF2", "BND" = "#80BDE9")) +
  theme(axis.line = element_line(color = 'black'),
        legend.position = "top",axis.text.y = element_text(color = 'black', size = 12),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())


# Extended Fig.1jklm: LRS TR mutation pattern ----------------------------------------

#input data include 3 column: TRID TR_Type Mean_TR_CNV_across_GTOP_samples
ggplot(LRS, aes(x = TR_Type, y = Mean_TR_Number,fill=TR_Type)) +
  geom_boxplot2(width = .7, width.errorbar = .5) +  
  labs(x = "Repeat units Length(bp)", y = "Mean number of RU copies at loci") +
  theme_classic() +
  theme(legend.position = "none")+
  scale_fill_manual(values = ru_color)+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(color = "black", size = 12), 
        axis.ticks = element_line(color = "black"))

#input data include 3 column: TRID TR_Type major_AF
ggplot(LRS,aes(major_AF)) +
  geom_histogram(binwidth = 0.01, fill = "#A2B5CD", color = "white")+
  scale_y_log10( breaks = c(1, 10, 100, 1000, 10000, 100000))+
  theme_classic()+
  labs(x = "Major allele frequency", y = "# of TR")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(color = "black", size = 12), 
        axis.ticks = element_line(color = "black"))

#input data include 3 column: TRID TR_Type allele_count
ggplot(LRS,aes(x=TR_Type,y=allele_unique_count,fill=TR_Type)) +
  geom_boxplot2(width = .7, width.errorbar = .5)+ 
  theme_classic()+
  scale_fill_manual(values = ru_color)+
  labs(y = "# of alleles per locus", x = "TR unit length")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(color = "black", size = 12), 
        axis.ticks = element_line(color = "black"),
        legend.position = "none")

#TR per genome count
LRS <- fread("Supp_Fig5k.TR_perGenome.txt")
LRS2 <- LRS %>%
  mutate(total = RU2 + RU3 + RU4 + RU5 + RU6 + VNTR) %>%
  arrange(desc(total)) %>%   
  mutate(Sample_ID = factor(Sample_ID, levels = Sample_ID)) %>%  
  select(Sample_ID, RU2, RU3, RU4,RU5,RU6, VNTR)

lrs_ind_count <- LRS2 %>%pivot_longer(cols = c("RU2","RU3","RU4","RU5","RU6","VNTR"),
               names_to = "TR_type", values_to = "count") %>% 
  mutate(TR_type=str_replace(TR_type,"RU",""))

ggplot(lrs_ind_count, aes(x = Sample_ID,   y = count, fill = TR_type)) +
  geom_col(position = "stack", width = 1) +  
  labs(x = "Sample", y = "# of TRs") +
  theme_classic() +
  theme(axis.line = element_line(color = 'black'),
        legend.position = "top",axis.text.y = element_text(color = 'black', size = 12),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())


