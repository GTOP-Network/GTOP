#===============================================#
# LRS read QC #
# Supp Figure 4   #
#===============================================#
library(ggplot2)
library(ggpubr)


setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig4/input")


# Supp.Figure.4a scatter plot of aligned reads and all reads --------------


df <- fread("supp4a.LR_read_QC_number.txt")

ggplot(df)+
  geom_point(aes(x=`Number of all reads`/1000000, y=`Number of aligned reads`/1000000,
                 color=Tissue))+
  scale_color_manual(values=setNames(df$Tissue_Color_Code, df$Tissue))+
  xlab("Number of all reads (million)")+
  ylab("Number of aligned reads (million)")+
  theme_pubr()

# Supp.Figure.4b distribution of Qscore -----------------------------------


df <- fread("supp8b.LR_read_QC_Qscore.txt.gz")
plot_df <- df %>% 
  pivot_longer(cols = everything(), names_to = "sample", values_to = "Qscore") %>% 
  separate(sample, into = c("sample", "tissuecolor"), sep = ";")
color <- plot_df %>% distinct(tissuecolor)

ggplot(plot_df)+
  geom_density(aes(x=Qscore, color=tissuecolor))+
  scale_color_manual(values=setNames(color$tissuecolor, color$tissuecolor))+
  xlab("Phred quality score")+
  ylab("Density")+
  theme_pubr()+
  theme(legend.position = "none")


# Supp.Figure.4c read length per tissue -----------------------------------


df <- fread("supp4c.LR_read_QC_tissue_read_length.txt.gz")
plot_df <- df %>% 
  pivot_longer(cols = everything(), names_to = "tissue", values_to = "readlength") %>% 
  separate(tissue, into = c("tissue", "tissuecolor"), sep = ";")
color <- distinct(plot_df, tissue, tissuecolor)

ggplot(plot_df)+
  geom_boxplot(aes(x=tissue, y=log10(readlength), fill=tissue), outlier.color = "NA")+
  scale_fill_manual(values=setNames(color$tissuecolor, color$tissue))+
  scale_y_continuous(limits = c(0, 4))+
  xlab("")+
  ylab("log10(read length)")+
  theme_pubr()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none")
