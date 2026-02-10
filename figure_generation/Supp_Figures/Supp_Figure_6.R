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
library(ggrastr)
library(ggupset)

setwd("/path/to/GTOP_code/supp/supp_fig6/input")
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

p1
#Supp.Fig.6b:  SV，
ggplot(result[result$type!="SNP",], aes(x = type, y = count, fill = category)) + 
  geom_bar(stat="identity", color = "black", width = 1) + 
  theme_classic() + 
  theme(axis.text = element_text(color = "black", size = 12), 
        axis.ticks = element_line(color = "black"),
        legend.position = "top") + 
  facet_wrap(~coverage, ncol = length(unique(result$coverage))) 

#Supp.Fig.6c:  TR
number <- fread("Supp_Fig6c.TR_downsample.txt",header = T)
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

# Supp.Fig.6d: small variant genotype accuarcy between LRS and SRS  --------------------------------------------------------
#The input file contains allele frequencies (AF) for overlapping variant sites between SRS and LRS. 
#The file is too large.
smoothScatter(dat$ALT_AF_LRS, dat$ALT_AF_SRS,
              nrpoints = 0, 
              bandwidth = 0.01,  
              xlab = "Allele Frequency in GTOP LRS",
              ylab = "Allele Frequency in GTOP SRS",
              main = "")

# Supp.Fig.6e: LRS specific small variant in difficult regions  --------------------------------------------------------
data <- tribble(~tech, ~type, ~variant, ~intersect, ~total,
                "LRS","specific","smallVariant",1642932+1770270,1806159+1789021,
                "LRS","share","smallVariant",2520341+996826,14485260+1690352,
                "SRS","specific","smallVariant",280861+1128651,697419+1313748,
                "SRS","share","smallVariant",2520341+996826,14485260+1690352)%>%
  mutate(proportion = intersect / total)%>% mutate(type=factor(type,level=c("specific","share")))%>%
  select(tech, type, variant, intersect, total) %>%
  pivot_wider(names_from = type, values_from = c(intersect, total), names_sep = "_") %>%
  mutate(specific_in_intersect = intersect_specific,
         specific_remaining = total_specific - intersect_specific,
         share_value = total_share) %>%
  select(tech, variant, specific_in_intersect, specific_remaining, share_value) %>%
  pivot_longer(cols = c(specific_in_intersect, specific_remaining, share_value),
               names_to = "category",
               values_to = "value") %>%
  group_by(tech,variant) %>%
  mutate(prop = value / sum(value))%>%
  mutate(category=factor(category,level=c("specific_in_intersect","specific_remaining","share_value")))

ggplot(data, aes(x = tech, y = prop, fill = category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(
    "specific_in_intersect" = "#9a3929",
    "specific_remaining" = "grey",
    "share_value" = "#8196b7" )) +
  ylab("Proportion") +
  xlab("Tech") +
  theme_classic()

# Supp.Fig.6f: SV number for deepvariant vs clair3  --------------------------------------------------------
library(VennDiagram)
draw.pairwise.venn(
  area1 = 16292086,     # total in DeepVariant
  area2 = 21053049,      # total in Clair3
  cross.area = 15959765, # overlap
  category = c("DeepVariant", "Clair3"),
  fill = c("#377EB8","#cdc2c9"),
  alpha = 0.6,
  lwd = 2,
  cex = 1.5,
  cat.cex = 1.4
)
# Supp.Fig.6g: SV compare in sniffles2/pbsv/cuteSV  --------------------------------------------------------

upset_data <- data.frame(
  combination = c(
    "Sniffles,PBSV,cuteSV",
    "Sniffles,PBSV",
    "Sniffles,cuteSV",
    "PBSV,cuteSV",
    "Sniffles",
    "PBSV",
    "cuteSV" ),
  count = c(72131, 11728, 8123, 6571, 8821, 3859, 4709)) %>%
  mutate(
    tools = str_split(combination, ","),
    label = case_when(
      combination == "Sniffles,PBSV,cuteSV" ~ "All Three",
      combination == "Sniffles,PBSV" ~ "Sniffles + PBSV",
      combination == "Sniffles,cuteSV" ~ "Sniffles + cuteSV",
      combination == "PBSV,cuteSV" ~ "PBSV + cuteSV",
      combination == "Sniffles" ~ "Sniffles only",
      combination == "PBSV" ~ "PBSV only",
      combination == "cuteSV" ~ "cuteSV only")) %>%
  arrange(desc(count)) %>%
  mutate(label = factor(label, levels = label),
         tools = str_split(combination, ",")  )


upset_data <- upset_data %>%
  mutate( display_label = case_when(
    combination == "Sniffles,PBSV,cuteSV" ~ "All Three",
    combination == "Sniffles,PBSV" ~ "Sniffles + PBSV",
    combination == "Sniffles,cuteSV" ~ "Sniffles + cuteSV",
    combination == "PBSV,cuteSV" ~ "PBSV + cuteSV",
    combination == "Sniffles" ~ "Sniffles only",
    combination == "PBSV" ~ "PBSV only",
    combination == "cuteSV" ~ "cuteSV only" )) %>%
  arrange(desc(count)) %>%
  mutate(display_label = factor(display_label, levels = display_label),
         sniffles = ifelse(str_detect(combination, "Sniffles"), 1, 0),
         pbsv = ifelse(str_detect(combination, "PBSV"), 1, 0),
         cutesv = ifelse(str_detect(combination, "cuteSV"), 1, 0))

point_data <- upset_data %>%
  pivot_longer(cols = c(sniffles, pbsv, cutesv), names_to = "tool", values_to = "present") %>%
  filter(present == 1) %>%
  mutate(tool = factor(tool,levels = c("sniffles", "pbsv", "cutesv"),
                       labels = c("Sniffles", "PBSV", "cuteSV")),
         tool_y = case_when( tool == "Sniffles" ~ 1,tool == "PBSV" ~ 2, tool == "cuteSV" ~ 3 ) )

p_points <- ggplot(point_data,
                   aes(x = display_label, y = tool_y)) +
  geom_point(size = 6, shape = 15, color = "black", fill = "black") +
  scale_y_continuous(
    breaks = 1:3,
    labels = c("Sniffles", "PBSV", "cuteSV"),
    limits = c(0.5, 3.5)) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 0, l = 10))

p_bars <- ggplot(upset_data, 
                 aes(x = display_label, y = count)) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.7) +
  geom_text(aes(label = format(count, big.mark = ",")),
            vjust = -0.5, size = 4, fontface = "bold") +
  labs(x = "Tool Combination", y = "Number of Variants") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.title.y = element_text(size = 14),
    plot.margin = margin(t = 0, r = 10, b = 10, l = 10)) +
  scale_y_continuous(
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.15)),
    limits = c(0, max(upset_data$count) * 1.15))

tool_totals <- data.frame(
  tool = c("Sniffles", "PBSV", "cuteSV"),
  total = c(100803, 94289, 91534),
  color = c("#FF6B6B", "#4ECDC4", "#45B7D1"))
p_tools <- ggplot(tool_totals, aes(x = reorder(tool,-total), y = total, fill = tool)) +
  geom_bar(stat = "identity", width = 0.7)

p_combined <- cowplot::plot_grid(p_bars,p_points,p_tools,ncol = 1)

print(p_combined)


# Supp.Fig.6h: SV compare with reported datasets  --------------------------------------------------------

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


# Supp.Fig.6i-k: LRS SV mutation pattern ---------------------------------------------------

dat <- fread("Supp_Fig6i-j.LRS_SV_info.txt",header = T)
dat <- dat %>% setnames(c("chr","start","ID","length","maf","hwe","end")) %>%
  mutate(svtype=word(ID,start=3,end=3,sep="\\_"))

#Supp.Fig.6i: maf
ggplot(dat, aes(maf)) +
  geom_histogram(binwidth = 0.01, fill = "#A2B5CD", color = "black",alpha=0.8) + 
  labs(x = "Minor allele frequency (log scale)", y = "# of SVs") +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 12), 
        axis.ticks = element_line(color = "black"))

#Supp.Fig.6j: length
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

#Supp.Fig.6k: pergenome count

datt<-fread("Supp_Fig6k.txt")
datt <- datt %>% mutate(sampleID=factor(sampleID,
                                        level=c(datt %>%group_by(sampleID) %>%summarise(total_sv = sum(sv_type_number, na.rm = TRUE)) %>%
                                         arrange(-total_sv) %>%pull(sampleID))))

ggplot(datt, aes(x = sampleID,  y = sv_type_number,fill=type)) +
  geom_col(position = "stack", width = 1) +  
  labs(x = "Sample", y = "# of SVs") +
  theme_classic() +
  scale_fill_manual(values = c("DEL" = '#5C87A6', "INS" = '#C56364', "DUP" = '#EEA3A3',
                               "INV" = "#BDDEF2", "BND" = "#80BDE9")) +
  theme(axis.line = element_line(color = 'black'),
        legend.position = "top",axis.text.y = element_text(color = 'black', size = 12),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())


# Supp.Fig.6l-n lmno: LRS TR mutation pattern ----------------------------------------

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

# Supp.Fig.6o -------------------------------------------------------------

#TR per genome count
LRS <- fread("Supp_Fig6o.TR_perGenome.txt")
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


