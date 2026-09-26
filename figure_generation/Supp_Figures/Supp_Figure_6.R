#==============================================#
# LRS & SRS Small variant/SV/TR#
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

source("geom_boxplot2.R")
setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig2567/input")


# Supp.Fig.6a: Small variant number for deepvariant vs clair3  --------------------------------------------------------
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

# Supp.Fig.6b-c: LRS specific small variant in difficult regions  --------------------------------------------------------

dat <- fread("./Supp_Fig6b-c.LRS_SRS_snv_compare.txt")
dat <- dat[grepl("difficult|Genome",dat$Region),]%>%
  mutate(Category=paste0(VarType_new,":",ifelse(grepl("_alldifficult",Region),"difficult",
                                                ifelse(grepl("_notinalldifficultregions",Region),"easy",
                                                       "Genome-wide"))))%>%
  dplyr::filter(!(VarType_new) %in% "INDEL_all")

make_plot <- function(dat, cols, name, keep){
  d <- melt(dat, id.vars=c("Category"), measure.vars=keep, variable.name="Class", value.name="Count")
  d$Class <- factor(d$Class, levels=keep)
  d$Prop <- d$Count / ave(d$Count, d$Category, FUN=sum)
  d$Category <- factor(d$Category,
                       levels = c("SNV:Genome-wide","INDEL_1_5:Genome-wide","INDEL_6_15:Genome-wide","INDEL_16:Genome-wide",
                                  "SNV:easy","INDEL_1_5:easy","INDEL_6_15:easy","INDEL_16:easy",
                                  "SNV:difficult","INDEL_1_5:difficult","INDEL_6_15:difficult","INDEL_16:difficult"))
  p1 <- ggplot(d, aes(Category, Count, fill=Class)) + geom_col(width=0.8) +
    scale_fill_manual(values=cols) + theme_classic(base_size=13) +
    labs(x=NULL, y=paste0(name," count")) +
    scale_y_continuous(limits = c(0,16500000))+ 
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="top")
  p2 <- ggplot(d, aes(Category, Prop, fill=Class)) + geom_col(width=0.8) +
    scale_fill_manual(values=cols) + scale_y_continuous(expand=c(0,0)) +
    theme_classic(base_size=13) +
    labs(x=NULL, y=paste0(name," proportion")) +
    theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="top")
  return(p1 / p2)
}

cols_lrs <- c(shared="#c9c9cb", LRS_specific="#8c1c21")
cols_srs <- c(shared="#c9c9cb", SRS_specific="#4e7266")

p_lrs <- make_plot(dat, cols_lrs, "LRS", c("LRS_specific","shared"));p_lrs
p_srs <- make_plot(dat, cols_srs, "SRS", c("SRS_specific","shared"));p_srs


# Supp.Fig.6d: small variant genotype accuarcy between LRS and SRS  --------------------------------------------------------
#The input file contains allele frequencies (AF) for overlapping variant sites between SRS and LRS. 
smoothScatter(dat$ALT_AF_LRS, dat$ALT_AF_SRS,
              nrpoints = 0, 
              bandwidth = 0.01,  
              xlab = "Allele Frequency in GTOP LRS",
              ylab = "Allele Frequency in GTOP SRS",
              main = "")


# Supp.Fig.6e: SV compare in sniffles2/pbsv/cuteSV  --------------------------------------------------------

upset_data <- fread("Supp_Fig6e.LRS_SV_tools_compare.txt")
upset_data <- upset_data %>%
  mutate( combination=factor(
    combination,
    levels=c("sniffles,pbsv,cutesv","pbsv","cutesv", "sniffles,cutesv","pbsv,cutesv","sniffles","sniffles,pbsv")),
    sniffles=ifelse(str_detect(combination,"sniffles"),1,0),
    pbsv=ifelse(str_detect(combination,"pbsv"),1,0),
    cutesv=ifelse(str_detect(combination,"cutesv"),1,0) )

point_data <- upset_data %>%
  pivot_longer(cols=c(sniffles,pbsv,cutesv),names_to="tool",values_to="present" ) %>%
  filter(present==1) %>%
  mutate(tool=factor(tool,levels=c("sniffles","pbsv","cutesv"),labels=c("sniffles","pbsv","cutesv")),y=as.numeric(tool))

p1 <- ggplot(point_data,aes(x=combination,y=y))+
  geom_point(size=5)+
  scale_y_continuous(breaks=1:3,labels=c("sniffles","pbsv","cutesv") )+
  theme_bw()+
  theme(axis.text.x=element_blank(),axis.title=element_blank(),panel.grid=element_blank())

p2 <- ggplot(upset_data,aes(x=combination,y=count))+
  geom_col(fill="steelblue")+
  scale_y_continuous(labels=scales::comma,expand=expansion(mult=c(0,0.15)))+
  theme_classic()+
  theme(axis.text.x=element_text(angle=45,hjust=1) )

count <- data.frame(tool=c("sniffles","pbsv","cutesv"),count=c( 102509, 95717, 91564))%>%
  mutate(tool=factor(tool,level=c("sniffles","pbsv","cutesv")))
p3 <- ggplot(count,aes(x=tool,y=count))+
  geom_col(fill="grey40")+
  theme_classic()

cowplot::plot_grid(p2,p1,p3,ncol=1,rel_heights=c(2,1,1))


# Supp.Fig.6f: LRS-SRS shared SV number in different BP distance --------------------------------------------------------

dat <- fread("Supp_Fig6f.LRS_SRS_SV_saturation.txt")
ggplot(dat,aes(x = threshold,y = 1-sensitivity, color = svtype, shape = svtype)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  theme_classic() +
  scale_color_manual(values = c("INS" = "#913628", "DEL" = "#227e85"  )) +
  scale_shape_manual(values = c( "INS" = 16,  "DEL" = 17)) +
  scale_x_continuous(limits = c(0, 500),breaks = seq(0, 500, by = 100)) +
  labs( x = "Breakpoint distance threshold (bp)",
        y = "Recall",color = "Dataset", shape = "SV type") +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.line = element_line(color = "black"),
        legend.position = "right")+
  scale_y_continuous(limits = c(0,1))

# Supp.Fig.6g: LRS/ SRS SV number ---------------------------------------------------

dat <- data.frame(tec = rep(c("LRS", "SRS"), each = 4),
                  TR_type = rep(c("INS", "DEL"), each = 2, times = 2),
                  group = rep(c("Shared", "Specific"), 4),
                  count = c(13097, 47293,11207, 14594,
                            8068, 1533,11226, 4546))

col_overlap  <- "#8da0cb";col_lrs_spec <- "#983628";col_srs_spec <- "#247E85"
fill_map <- c("Shared.LRS"   = col_overlap,"Shared.SRS"   = col_overlap,"Specific.LRS" = col_lrs_spec,"Specific.SRS" = col_srs_spec)
ggplot() +
  geom_col(data = dat %>% filter(tec == "LRS") %>%mutate(group=factor(group,level=c("Specific","Shared"))),
           aes(y = count, x = TR_type, fill = interaction(group, tec)), position = "stack") +
  geom_col(data = dat %>% filter(tec == "SRS") %>%mutate(group=factor(group,level=c("Specific","Shared"))),
           aes(y = -count, x = TR_type, fill = interaction(group, tec)), position = "stack") +
  scale_fill_manual(values = fill_map, name = "Group") +
  scale_y_continuous(labels = function(x) {
    lx <- abs(x)
    ifelse(lx >= 1000, paste0(lx / 1000, "k"), lx) },
    limits = c(-60390 * 1.1, 60390 * 1.1) ) + # max total number
  geom_hline(yintercept = 0, color = "black") +
  labs(x = "", y = "Number of SVs") +
  theme_classic(base_size = 13) +
  theme(axis.text = element_text(size = 11),
        legend.position = "top")



# Supp.Fig.6h LRS / SV length distribution ---------------------------------


ggplot(fread("Supp_Fig6h.LRS_SRS_SV_DEL_length.txt"),aes(x=Var1,y=Freq,color=Var2)) + geom_line(size=1.5) + theme_pubr() + 
  scale_y_log10(breaks=c(1,10,100,1000,10000),labels=c("1","10","100","1K","10K")) + 
  scale_x_log10(breaks=c(100,1000,10000,100000,500000,100000000),labels=c("100","1K","10K","100K","500K",">=1M")) + 
  scale_color_manual(breaks = c("LRS","SRS"),values = c("#9c3929","#1c9099"))

ggplot(fread("Supp_Fig6h.LRS_SRS_SV_INS_length.txt"),aes(x=Var1,y=Freq,color=Var2)) + geom_line(size=1.5) + theme_pubr() + 
  scale_y_log10(breaks=c(1,10,100,1000,10000),labels=c("1","10","100","1K","10K")) + 
  scale_x_log10(breaks=c(100,1000,10000,100000),labels=c("100","1K","10K","100K")) + 
  scale_color_manual(breaks = c("LRS","SRS"),values = c("#9c3929","#1c9099"))


# Supp.Fig.6i LRS specific SV enrichment ---------------------------------

dat <- fread("Supp_Fig6i.LRS_specific_SV_enrichment.txt")%>%
  dplyr::filter(tech %in% c("LRS"))  %>% 
  mutate(region=factor(region,level=c("tr","low_mapp","sd")))%>%
  mutate(group=factor(group,level=c("specific","shared")))

ggplot(dat,aes(x = region, y = prop, fill = group)) +
  geom_col(position = position_dodge(0.9), width = 0.8) +
  geom_text(data = subset(dat,group == "specific"),
            aes(x = region, y = prop + 0.03, label = sig_label),
            position = position_dodge(0.9),
            size = 5,vjust = 0) +
  labs( x = "", y = "Proportion of SVs within region") +
  scale_fill_manual(values = c("shared" = "#8da0cb", "specific" = "#9c3929"), name = "SV group") +
  theme_classic(base_size = 14) +
  theme( legend.position = "top",
         strip.text = element_text(face = "bold", size = 13),
         axis.text.x = element_text(angle = 0, hjust = 0.5))



