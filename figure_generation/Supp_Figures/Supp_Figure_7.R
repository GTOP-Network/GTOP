#==============================================#
# Small variant /SV /TR #
# Supp-Figure-7#
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
library(ggrastr)
library(ggupset)

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig2567/input")
source("geom_boxplot2.R")

# Supp.Fig.7a-c:  Small variant /SV /TR number in diff LRS coverage ------------------------------------------------------------------

#Supp.Fig.7a:  Small variant
number <- fread("Supp_Fig7a.small_variant_downsample.txt")
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

p1

# Supp.Fig.7b:  SV -------------------------------------------------------


dat <- fread("Supp_Fig7b.sv_downsample.txt")
dat$Depth <- factor(dat$Depth,levels = c("5","10","15","20","25","30","35","40"))
dat$SVTYPE <- factor(dat$SVTYPE,levels = c("INV","DUP","BND","DEL","INS"))
ggplot(dat, aes(x=Depth, y=Count, fill=SVTYPE)) +
  geom_bar(stat="identity", position="stack", width=0.7) +
  labs(x="Sequencing depth", y="Number of SVs", fill="SV type") +
  theme_classic(base_size=14)+
  scale_fill_manual(values = c("DEL" = '#5C87A6', "INS" = '#C56364', "DUP" = '#EEA3A3',"INV" = "#BDDEF2", "BND" = "#80BDE9"))


# Supp.Fig.7c:  TR --------------------------------------------------------


number <- fread("Supp_Fig7c.TR_downsample.txt",header = T)
number$V3 <- c("5x","10x","15x","20x","25x","30x","35x","40x")
number$V3 <- factor(number$V3,levels = c("5x","10x","15x","20x","25x","30x","35x","40x"))
data_long <- number %>%pivot_longer(cols = c(STR_Count, VNTR_Count),names_to = "Type",values_to = "Count")
ggplot(data_long, aes(x = V3, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("STR_Count" = "#3e669a", "VNTR_Count" = "#ba9fad"),
                    labels = c("STR", "VNTR")) +
  labs() +
  theme_classic() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.position = "top")


# Supp.Fig.7d-f: LRS SV mutation pattern ---------------------------------------------------

dat <- fread("Supp_Fig7d-f.LRS_SV_info.txt",header = T)

#Supp.Fig.7d: maf
ggplot(dat, aes(maf)) +
  geom_histogram(binwidth = 0.01, fill = "#A2B5CD", color = "black",alpha=0.8) + 
  labs(x = "Minor allele frequency (log scale)", y = "# of SVs") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 12), 
        axis.ticks = element_line(color = "black"))

#Supp.Fig.7e: length
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

#Supp.Fig.7f: pergenome count

dat<-fread("Supp_Fig7f.SV_PerSample_count.txt")
dat <- dat %>% mutate(sampleID=factor(sampleID,level=c(dat %>%group_by(sampleID) %>%
                                                         summarise(total_sv = sum(sv_type_number, na.rm = TRUE)) %>%
                                         arrange(-total_sv) %>%pull(sampleID))))

ggplot(dat, aes(x = sampleID,  y = sv_type_number,fill=type)) +
  geom_col(position = "stack", width = 1) +  
  labs(x = "Sample", y = "# of SVs") +
  theme_classic() +
  scale_fill_manual(values = c("DEL" = '#5C87A6', "INS" = '#C56364', "DUP" = '#EEA3A3',
                               "INV" = "#BDDEF2", "BND" = "#80BDE9")) +
  scale_y_continuous(limits = c(0, 25000), breaks = seq(0, 25000, by = 5000)) +
  theme(axis.line = element_line(color = 'black'),
        legend.position = "top",axis.text.y = element_text(color = 'black', size = 12),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Supp.Fig.7g -------------------------------------------------------------

ac_df <- fread("Supp_Fig7g.TR_allele_diversity.txt")
ggplot(ac_df, aes(x = factor(threshold), y = OR)) +
  geom_point(size = 3, color = "#993828") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "#993828") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x = "Allele count threshold",
    y = "Odds Ratio (LRS-specific vs Shared)") +
  theme_classic(base_size = 14)


# Supp.Fig.7h-k lmno: LRS TR mutation pattern ----------------------------------------

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

# Supp.Fig.7k -------------------------------------------------------------

#TR per genome count
LRS <- fread("Supp_Fig7k.TR_perGenome.txt")
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
  theme_classic() +scale_fill_igv()+
  theme(axis.line = element_line(color = 'black'),
        legend.position = "top",axis.text.y = element_text(color = 'black', size = 12),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())


