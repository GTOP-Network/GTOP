#====================================#
# Cross-tissue co-expression module#
# # Supp-Figure-12 # #
#===================================#

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig12/input")


# supp.Figure.12a proportion of novel & annotated in module ---------------

library(dplyr)
library(stringr)

load("all_gene_in_module.RData")
length(unique(all$gene_id))

all$type<-ifelse(str_detect(all$gene_id,"GTOP"),"novel","annotated")
table(all$type)

bar_plot<- all %>%
  dplyr::count(Tissue, type, name = "n")
tissue_order <- bar_plot %>%
  filter(type == "annotated") %>%
  arrange(desc(n)) %>%
  pull(Tissue)

bar_plot$Tissue <- factor(bar_plot$Tissue, levels = tissue_order)
bar_plot$n<-bar_plot$n/1000
library(ggplot2)
ggplot(bar_plot, aes(
  x = Tissue,
  y = n,
  fill = type
)) +
  scale_y_continuous(limits = c(0, 130),breaks = c(0,30,60,90,130)) +  
  geom_col(
    position = position_dodge(width = 0.7),
    width = 0.6
  ) +
  labs(
    x = NULL,
    y = "Number of transcripts(10³)",
    fill = NULL
  ) +
  scale_fill_manual(values = c("novel" = "#9d3929", "annotated" = "#7d8bad")) +  # 指定颜色
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 1
    )
  )


# Supp.Figure.12b annotated modules & unannotated modules -----------------

load("need_module.RData")
all_tissue_module<- all %>%
  distinct(Tissue, module) %>%   
  dplyr::count(Tissue, name = "n_modules")

all_go_module<-need_module %>%
  distinct(Tissue, Module) %>%  
  dplyr::count(Tissue, name = "n_modules")


plot_data<-all_tissue_module
plot_data$go_module<-all_go_module$n_modules
colnames(plot_data)<-c("Tissue","Unannotated_modules","Annotated_modules")

plot_data$Unannotated_modules<-plot_data$Unannotated_modules-plot_data$Annotated_modules

library(dplyr)
library(tidyr)

plot_long <- plot_data %>%
  mutate(Total = Annotated_modules + Unannotated_modules) %>%
  arrange(desc(Total)) %>%
  mutate(
    Tissue = factor(Tissue, levels = Tissue) 
  ) %>%
  pivot_longer(
    cols = c(Unannotated_modules, Annotated_modules),
    names_to = "Module_type",
    values_to = "Count"
  ) %>%
  mutate(
    Module_type = factor(
      Module_type,
      levels = c( "Annotated_modules","Unannotated_modules")
    )
  )
library(ggplot2)
library(ggbreak)
ggplot(plot_long, aes(x = Tissue, y = Count, fill = Module_type)) +
  geom_col(width = 0.8) +
  scale_fill_manual(
    values = c(
      "Unannotated_modules" = "#BBBBBB",  
      "Annotated_modules"   = "#7d8bad"
    )
  ) +
  labs(
    x = NULL,
    y = "Number of modules",
    fill = NULL
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.y.right = element_blank(),
    axis.text.y.right = element_blank()
  )+
  #scale_y_break(c(130, 230)) +   
  scale_y_continuous(
    breaks = c(0, 50, 100, 200,300, 400, 500)  
  )
# Supp.Figure.13c ---------------------------------------------------------

# group1:gene with 1 isofrom
# group2:gene with multiple isofrom all in same module
# group3:gene with multiple isofrom across multiple modules

load("all_results.RData")
plot_df<-unique(all_results[,3:5])
plot_df <- data.frame(table(plot_df$Tissue,plot_df$group))
colnames(plot_df)<-c("tissue","group","count")

library(dplyr)

plot_df <- plot_df %>%
  group_by(tissue) %>%
  mutate(total_count = sum(count)) %>%
  ungroup() %>%
  mutate(
    tissue = factor(tissue,
                    levels = unique(tissue[order(-total_count)]))
  )

plot_df$count<-plot_df$count/1000
ggplot(plot_df, aes(x = tissue, y = count, fill = group)) +
  geom_col(width = 0.8) +
  scale_fill_manual(
    values = c(
      "g1" = "#87a2cb",  
      "g2"   = "#8cae6a",
      "g3"   = "#c9c9c9"
    )
  ) +
  labs(
    x = NULL,
    y = "Number of genes (10³)",
    fill = NULL
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
    #axis.line.x = element_blank()
  )

# supp.Figure.13d WGCNA for heart novel ----------------------------------

load("heart_go_use_novel.RData")
load("merged_expr_novel.RData")
load("sample_annot_full_1586.RData")
sample_annot_full$sample_id<-gsub("AGTEX","GTOP",sample_annot_full$sample_id)
heart_sample<-sample_annot_full[str_detect(sample_annot_full$Tissue,"Heart"),]

tissue <- heart_sample$Tissue
names(tissue) <- heart_sample$sample_id

# Z-score
merged_expr_norm <- t(scale(t(merged_expr), center = TRUE, scale = TRUE))
merged_expr_norm <- as.data.frame(merged_expr_norm)

my_colors <- c(colorRampPalette(c("#8089a9", "white"))(50),
               colorRampPalette(c("white", "#bb4633"))(50))
my_breaks <- c(seq(-2, 0, length.out = 51), seq(0.01, 2, length.out = 50))

library(pheatmap)
p<-pheatmap(merged_expr_norm,
            color = my_colors,
            breaks = my_breaks,
            cluster_rows = T,
            cluster_cols = T,
            border_color = NA,
            show_rownames = T,
            show_colnames = T
)

p

row_order <- p$tree_row$order
row_names_sorted <- rownames(merged_expr_norm)[row_order]

go_result_bar$log10p<- -log10(go_result_bar$p.adjust)
library(ggplot2)
go_result_bar<-go_result_bar[match(row_names_sorted,go_result_bar$Module),]
go_result_bar$Description <- factor(go_result_bar$Description, levels = rev(go_result_bar$Description))

ggplot(go_result_bar, aes(x = Description, y = log10p)) +
  geom_bar(stat = "identity", width = 0.7, fill = "#56B4E9") +  
  coord_flip() +  
  labs(x = "GO Term", y = "-log10(p.adjust)", title = "GO Enrichment") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),          
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )




