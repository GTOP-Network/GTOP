#==============================================#
# Heritability and LRS enhance disease coloc #
# Supp-Figure-40#
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
library(UpSetR)

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig40.R4.Heritability-long-reads-coloc/")

# Supp.Fig.40a: heritability --------------------------------------------------------
df <- fread("./input/Supp_Fig40a.txt")
sumdf <- fread("./input/Supp_Fig40a.allsum.txt")

df$QTL<-factor(df$QTL,levels = c("eQTL","juQTL","tuQTL"))
df$Category<-factor(df$Category,levels = c("Musculoskeletal","Quantitative_trait","Dermatologic",
                                           "Digestive","Circulatory_system","Endocrine_metabolic",
                                           "Respiratory","Neoplasms","Genitourinary"))
ggplot(
  df ,aes(  x = Category,  y = median_Prop_h2,  fill = QTL)
) +
  geom_col(  position = position_dodge(width = 0.75),  width = 0.7) +
  geom_errorbar(  aes(    ymin = q1,    ymax = q3  ),  position = position_dodge(width = 0.75),  width = 0.2,  linewidth = 0.4) +
  geom_hline(data=sumdf, 
  aes(yintercept = median, color=QTL, fill=QTL),
  linetype = "dashed",linewidth = 0.5) +
  scale_fill_manual(values=c("eQTL"="#95c595","juQTL"="#5181b1","tuQTL"="#c97575"))+
  scale_color_manual(values=c("eQTL"="#95c595","juQTL"="#5181b1","tuQTL"="#c97575"))+
  theme_pubr() +
  labs(  x = "Trait category",  y = "Median Prop h²",  fill = "QTL type") +
  theme(panel.grid.major.x = element_blank(),axis.text.x = element_text(  angle = 45,  hjust = 1,  vjust = 1))


# Supp.Fig.40b-d: Disease Colocalized by Novel gene and transcript  --------------------------------------------------------

trait_colocrs=c(
  "Quantitative_trait" = "#2F4F4F",
  "Circulatory_system" = "#487c51",
  "Dermatologic" = "#7690a4",
  "Digestive" = "#f0d795",
  "Endocrine_metabolic" = "#4682B4",
  "Genitourinary" = "#B0C4DE",
  "Infectious_diseases" = "#9e4832",
  "Mental_disorders" = "#7c776b",
  "Musculoskeletal" = "#DB7093",
  "Neoplasms" = "#7371a3",
  "Neurological" = "#00008B",
  "Respiratory" = "#5192c1",
  "Sense_organs" = "#c4598a",
  "Symptoms" = "#2F4F4F",
  "Hematopoietic"="#00bfc4",
  "eQTL" = "#9ac294",
  "sQTL" = "#7b86a7",
  "eQTL,sQTL" = "#bc9cc1")

plot_xqtl_class <- function(x, name, class_colors){
  
  plot_df <- x %>%
    group_by(phecode_abbr, trait_class, locus_id) %>%
    summarise(Class_group = case_when(
      all(Class == "Known") ~ "Known_only",
      all(Class == "Novel") ~ "Novel_only",
      TRUE ~ "Both"), .groups="drop") %>%
    group_by(phecode_abbr, trait_class, Class_group) %>%
    summarise(count=n_distinct(locus_id), .groups="drop") %>%
    group_by(phecode_abbr, trait_class) %>%
    mutate(gwas_total=sum(count)) %>%
    ungroup() %>%
    group_by(trait_class) %>%
    mutate(trait_total=sum(gwas_total)) %>%
    ungroup()
  
  prop <- plot_df %>% group_by(Class_group) %>%
    summarise(n=sum(count), .groups="drop") %>%
    mutate(proportion=n/sum(n), xqtl=name)
  
  print(prop)
  
  x_order <- plot_df %>%
    select(phecode_abbr, trait_class, gwas_total, trait_total) %>%
    distinct() %>%
    arrange(desc(trait_total), desc(gwas_total)) %>%
    pull(phecode_abbr)
  
  plot_df <- plot_df %>%left_join(data.frame(phecode_abbr=x_order,pos=seq_along(x_order)),by="phecode_abbr")
  plot_df$Class_group <- factor(plot_df$Class_group,levels=c("Known_only","Both","Novel_only"))
  
  p <- ggplot() +
    geom_col(data=plot_df,aes(x=pos,y=count,fill=Class_group),width=0.8,position="stack") +
    scale_fill_manual(values=class_colors,name="Class") +
    ggnewscale::new_scale_fill() +
    geom_rect(data=plot_df %>%group_by(trait_class) %>%
                summarise(xmin=min(pos)-0.4,xmax=max(pos)+0.4,.groups="drop"),
              aes(xmin=xmin,xmax=xmax,ymin=-max(plot_df$count,na.rm=T)*0.12,
                  ymax=-max(plot_df$count,na.rm=T)*0.06,fill=trait_class)) +
    scale_fill_manual(values=trait_colocrs,name="Category") +
    scale_x_continuous(expand=c(0,0),breaks=NULL,labels=NULL) +
    scale_y_continuous(
      expand=expansion(mult=c(0.05,0.05)),
      limits=c(-max(plot_df$count,na.rm=T)*0.12,
               max(plot_df$gwas_total,na.rm=T)*1.05)) +
    labs(title=name, x=NULL,y="Number of colocalized GWAS loci") +
    theme_classic() +
    theme(axis.text.x=element_blank(),legend.position="right")
  
  return(list(plot=p, proportion=prop, data=plot_df))
}

eqtl_color <-  c("Known_only"="#ddecdb","Both"="#00ba38","Novel_only"="#487c51")
juqtl_color <- c("Known_only"="#c1c8d8","Both"="#4682b4","Novel_only"="#0f4a75")
tuqtl_color <- c("Known_only"="#eed9d8","Both"="#db7093","Novel_only"="#9d3929")


## eQTL
eqtl <- fread("./input/Supp_Fig40.novel_eqtl_coloc.txt")
plot_xqtl_class(eqtl,"eQTL",eqtl_color)

# junction usage sQTL
juqtl <- fread("./input/Supp_Fig40.novel_juqtl_coloc.txt")
plot_xqtl_class(juqtl,"juQTL",juqtl_color)

# transcript usage sQTL
tuqtl <- fread("./input/Supp_Fig40.novel_tuqtl_coloc.txt")
plot_xqtl_class(tuqtl,"tuQTL",tuqtl_color)

# Supp.Fig.40e-gtissue specific coloc ---------------------------------------------------
tissue_color <- read.csv("input/GTOP_tissue_coloc_code",header = T)%>%
  mutate(Tissue_Color_Code=paste0("#",Tissue_Color_Code))
color_vec <- setNames(tissue_color$Tissue_Color_Code, tissue_color$Tissue)

## eQTL
e_upset_matrix <- fread("./input/Supp_Fig40.eqtl_tissue_coloc.txt")
bar_colors <- sapply(colnames(e_upset_matrix), function(t) ifelse(t %in% names(color_vec), color_vec[t], "#CCCCCC"))
eqtl_tissues_ordered <- names(sort(colSums(e_upset_matrix), decreasing = T))
eqtl_bar_colors_ordered <- bar_colors[eqtl_tissues_ordered]
upset(e_upset_matrix,
      sets = eqtl_tissues_ordered,
      sets.bar.color = eqtl_bar_colors_ordered,
      order.by = "freq",
      empty.intersections = "on",
      main.bar.color = "grey",mb.ratio = c(0.4,0.6),
      show.numbers = FALSE)

# junction usage sQTL
ju_upset_matrix <- fread("input/Supp_Fig40.juqtl_tissue_coloc.txt")
bar_colors <- sapply(colnames(ju_upset_matrix), function(t) ifelse(t %in% names(color_vec), color_vec[t], "#CCCCCC"))
juqtl_tissues_ordered <- names(sort(colSums(ju_upset_matrix), decreasing = T))
juqtl_bar_colors_ordered <- bar_colors[juqtl_tissues_ordered]
upset(ju_upset_matrix,
      sets = juqtl_tissues_ordered,
      sets.bar.color = juqtl_bar_colors_ordered,
      order.by = "freq",
      empty.intersections = "on",
      main.bar.color = "grey",mb.ratio = c(0.4,0.6),
      show.numbers = FALSE)

# transcript usage sQTL
tu_upset_matrix <- fread("input/Supp_Fig40.tuqtl_tissue_coloc.txt")
bar_colors <- sapply(colnames(tu_upset_matrix), function(t) ifelse(t %in% names(color_vec), color_vec[t], "#CCCCCC"))
tuqtl_tissues_ordered <- names(sort(colSums(tu_upset_matrix), decreasing = T))
tuqtl_bar_colors_ordered <- bar_colors[tuqtl_tissues_ordered]
upset(tu_upset_matrix,
      sets = tuqtl_tissues_ordered,
      sets.bar.color = tuqtl_bar_colors_ordered,
      order.by = "freq",
      empty.intersections = "on",
      main.bar.color = "grey",mb.ratio = c(0.4,0.6),
      show.numbers = FALSE)



