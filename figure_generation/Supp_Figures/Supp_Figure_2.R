#==============================================#
# LRS WGS QC #
# Supp-Figure-2#
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
library(ggExtra)

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig2567/input")


# Supp.Fig.2a-c:  LRS WGS QC Statistics ------------------------------------------------------------------

dat <- fread("Supp_Fig2.LRS_QC.txt")
dat$HiFi_Read_Length_N50 <- dat$HiFi_Read_Length_N50/1000
dat$base <- dat$HiFi_Yield/10^9
dat$coverage=dat$HiFi_Yield/3.1e9
dat$median_percent_identity <- (dat$Mapped_Reads/dat$HiFi_Reads)*100
p1 <- ggplot(dat, aes(x = coverage, y = HiFi_Read_Length_N50)) +
  geom_point(color = "black",size=3) +  
  geom_vline(xintercept =median(dat$coverage) , linetype = "dashed", color = "grey") + 
  geom_hline(yintercept = median(dat$HiFi_Read_Length_N50), linetype = "dashed", color = "grey") +  
  labs(x = "Sequencing Coverage(x)",y = "Read Length N50(Kb)") +
  theme_classic()
p_with_marginals3 <- ggMarginal(
  p1,type = "density",      
  fill = "skyblue",         
  alpha = 0.5,margins = "both")
print(p_with_marginals3)

p2 <- ggplot(dat, aes(x = Median_Q, y = median_percent_identity)) +
  geom_point(color = "black",size=3) +  
  geom_vline(xintercept =median(dat$Median_Q) , linetype = "dashed", color = "grey") + 
  geom_hline(yintercept = median(dat$median_percent_identity), linetype = "dashed", color = "grey") +  
  labs(x = "Median Read Qscore",y = "Mapping rate(%)") +
  theme_classic()+
  scale_y_continuous(limits = c(99,100))+
  scale_x_continuous(limits = c(30,35.5))
p_with_marginals3 <- ggMarginal(
  p2,type = "density",      
  fill = "pink",         
  alpha = 0.5,margins = "both")
print(p_with_marginals3)

p3 <- ggplot(dat,aes(x = GC)) +
  geom_histogram(binwidth = 0.05, fill = "grey", color = "black") +
  theme_classic() +
  geom_vline(xintercept = median(dat$GC), linetype = "dashed", size = 0.8, color ="#3d5488") + 
  annotate("text", x = median(dat$GC), y = max(hist(dat$GC, plot = FALSE)$counts)*0.1,
           label = paste0("Median: ", round(median(dat$GC), 2)), color = "black", fontface = "bold") +  
  labs(x = "GC content(%)",y = "Number of samples")

# Supp.Fig.2d-k: SRS WGS QC Statistics  --------------------------------------------------------
dat <- fread("Supp_Fig2.SRS_QC.txt")
coverage_plot <- ggplot(dat, aes(x = Coverage)) +
  geom_histogram(binwidth = 2, fill = "grey", color = "black") +
  theme_classic() +
  geom_vline(xintercept = median(dat$Coverage), linetype = "dashed", size = 0.8, color ="#3d5488") + 
  annotate("text", x = median(dat$Coverage), y = max(hist(dat$Coverage, plot = FALSE)$counts)*0.5,
           label = paste0("Median: ", round(median(dat$Coverage), 2)), color = "black", fontface = "bold") +  
  labs(x = "Coverage (X)",y = "Number of samples")+
  theme(axis.line=element_line(color='black'),
        axis.text=element_text(color='black'),
        legend.key=element_blank())

dat$total_bases_GB <- dat$after_filtering__total_bases / 1e9
q30_plot <- ggplot(dat, aes(x = total_bases_GB,y = after_filtering__q30_rate)) +
  geom_point(size = 2, color = "grey") +
  scale_y_continuous(limits = c(0.90,1))+
  geom_vline(xintercept = median(dat$total_bases_GB), linetype = "dashed", size = 0.8, color ="#3d5488") + 
  geom_hline(yintercept = median(dat$after_filtering__q30_rate), linetype = "dashed", size = 0.8, color ="#3d5488") + 
  theme_classic() +
  labs(x = "Sequencing(GB)",y = "Proportion of Q30")+
  theme(axis.line=element_line(color='black'),
        axis.text=element_text(color='black'),
        legend.key=element_blank())
median(dat$total_bases_GB)#145.31
median(dat$after_filtering__q30_rate)#0.967432

insert_plot <- ggplot(dat, aes(x = MEDIAN_INSERT_SIZE)) +
  geom_histogram(binwidth = 5, fill = "grey", color = "black") +
  theme_classic() +
  geom_vline(xintercept = median(dat$MEDIAN_INSERT_SIZE, na.rm = TRUE), linetype="dashed", color="#3d5488", size=0.8) +
  geom_vline(xintercept = median(dat$MEDIAN_INSERT_SIZE, na.rm = TRUE) + 4*sd(dat$MEDIAN_INSERT_SIZE, na.rm = TRUE),
             linetype="dashed", color="#3d5488", size=0.8) +
  geom_vline(xintercept = median(dat$MEDIAN_INSERT_SIZE, na.rm = TRUE) - 4*sd(dat$MEDIAN_INSERT_SIZE, na.rm = TRUE), linetype="dashed", color="#3d5488", size=0.8) +
  labs(x = "Median insert size (bp)",y = "Number of samples")+
  theme(axis.line=element_line(color='black'),
        axis.text=element_text(color='black'),
        legend.key=element_blank())

contamination_plot <- ggplot(dat, aes(x = `FREEMIX(Alpha)`)) +
  geom_histogram(binwidth = 0.0001, fill = "grey", color = "black") +
  theme_classic() +scale_x_continuous(limits = c(0,0.005))+
  labs(x = "Contamination rate",y = "Number of samples")+
  theme(axis.line=element_line(color='black'),
        axis.text=element_text(color='black'),
        legend.key=element_blank())

chimeric_plot <- ggplot(dat, aes(x = PCT_CHIMERAS)) +
  geom_histogram(binwidth = 0.005, fill = "grey", color = "black") +
  theme_classic() +scale_x_continuous(limits = c(0,0.08))+
  labs(x = "Percentage of chimeric reads",y = "Number of samples")+
  theme(axis.line=element_line(color='black'),
        axis.text=element_text(color='black'),
        legend.key=element_blank())

sex_plot <- ggplot(dat, aes(x = chrX_Autosome,y = chrY_Autosome)) +
  geom_point(size = 2,color = "grey") +
  theme_classic() +
  geom_vline(xintercept = 1, linetype="dashed", color="#3d5488", size=0.8) +
  geom_hline(yintercept = 0.5, linetype="dashed", color="#3d5488", size=0.8) +
  labs(x = "X/A depth ratio",y = "Y/A depth ratio")+
  theme(axis.line=element_line(color='black'),
        axis.text=element_text(color='black'),
        legend.key=element_blank())


mapping_plot <- ggplot(dat, aes(x = primary_proportion,y = properly_proportion)) +
  geom_point(size = 2,color = "grey") +
  theme_classic() +
  scale_x_continuous(limits = c(95,100))+
  scale_y_continuous(limits = c(90,100))+
  geom_vline(xintercept = median(dat$primary_proportion), linetype="dashed", color="#3d5488", size=0.8) +
  geom_hline(yintercept = median(dat$properly_proportion), linetype="dashed", color="#3d5488", size=0.8) +
  labs(x = "Primary mapping rate",y = "Properly mapping rate")+
  theme(axis.line=element_line(color='black'),
        axis.text=element_text(color='black'),
        legend.key=element_blank())

gc_plot <- ggplot(dat, aes(x = after_filtering__gc_content)) +
  geom_histogram(binwidth = 0.005, fill = "grey", color = "black") +
  theme_classic() +
  geom_vline(xintercept = median(dat$after_filtering__gc_content), linetype = "dashed", size = 0.8, color ="#3d5488") + 
  annotate("text", x = median(dat$after_filtering__gc_content), y = max(hist(dat$after_filtering__gc_content, plot = FALSE)$counts)*0.5,
           label = paste0("Median: ", round(median(dat$after_filtering__gc_content), 2)), color = "black", fontface = "bold") +  
  labs(x = "GC content(%)",y = "Number of samples")+
  theme(axis.line=element_line(color='black'),
        axis.text=element_text(color='black'),
        legend.key=element_blank())

# Supp.Fig.2l: Kinship and IBD ----------------------------------------------------------------------

dat <- fread("Supp_Fig2l.kinship_IBD.txt")
y0 <- min(dat$Kinship)
y1 <- max(dat$Kinship)
if (y1 < 0.1) {y1 <- 0.1}
p1 <- ggplot(dat, aes(x = Z0, y = Kinship)) +
  geom_point(color = "grey", size = 2) +
  scale_x_continuous(name = "Probability of zero IBD") +
  scale_y_continuous(
    name = "Kinship coefficient",
    limits = c(y0,y1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  theme_classic() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(color = "black", size = 10))






























