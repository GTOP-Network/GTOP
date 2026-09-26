#====================================#
# Cross-tissue co-expression module#
# # Supp-Figure-17 # #
#===================================#

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig17.figure_12/input")


# supp.Figure.17a  -----number of transcripts in module----------

library(dplyr)
library(stringr)
library(ggplot2)

load("suppfig.a.RData")

tissue_order <- plot_df_long %>%
  filter(category == "annotated") %>%
  arrange(desc(count)) %>%
  pull(Tissue)

plot_df_long$Tissue <- factor(plot_df_long$Tissue, levels = tissue_order)

# 如果count是原始数值，需要转换成千为单位
plot_df_long <- plot_df_long %>%
  mutate(count_k = count / 1000)

p_gtop <- ggplot(plot_df_long, aes(x = Tissue, y = count_k, fill = category)) +
  geom_col(width = 0.7, position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("annotated" = "#a1a7bd", "novel" = "#a0523e"),
                    labels = c("annotated" = "Annotated", "novel" = "Novel")) +
  labs(x = NULL, y = expression("Number of transcripts in modules (10"^3*")"), fill = NULL) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.85, 0.9),
    legend.background = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 0.5)
  )

p_gtop


# Supp.Figure.17b  --------number of modules------------------------------------------------

bar_width <- 0.5
offset <- 0
load("suppfig.b.RData")
tissue_order <- plot_long %>%
  group_by(Tissue) %>%
  summarise(total = sum(Count)) %>%
  arrange(desc(total)) %>%
  pull(Tissue)

plot_long$Tissue <- factor(plot_long$Tissue, levels = tissue_order)
plot_long$Fill_group<-factor(plot_long$Fill_group,levels = c("GTOP_Annotated_modules","GTOP_Unannotated_modules"))
n_module <- ggplot(plot_long, aes(x = Tissue, y = Count, fill = Fill_group)) +
  geom_col(width = 0.5) +
  scale_fill_manual(
    values = c(
      "GTOP_Unannotated_modules" = "#bfbfbf",
      "GTOP_Annotated_modules"   = "#7d8bad"
    ),
    breaks = c("GTOP_Annotated_modules", "GTOP_Unannotated_modules"),
    labels = c("Annotated modules", "Unannotated modules")
  ) +
  labs(x = NULL, y = "Number of modules", fill = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

n_module
# Supp.Figure.17c ---------------------------------------------------------
load("suppfig.c.RData")
ggplot(plot_df, aes(x = tissue, y = count, fill = group)) +
  geom_col(width = 0.8) +
  scale_fill_manual(
    values = c(
      "g1" = "#87a2cb",  
      "g2"   = "#8cae6a",
      "g3"   = "#c8c8c8"
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
# supp.Figure.17d  -------heatmap of novel transcripts in heart---------------------------

load("heart_go_use_novel.RData")
load("merged_expr_novel.RData")

my_colors <- c(colorRampPalette(c("#8089a9", "white"))(50),
               colorRampPalette(c("white", "#bb4633"))(50))
my_breaks <- c(seq(-1.5, 0, length.out = 51), seq(0.01, 1.5, length.out = 50))

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

p<-ggplot(go_result_bar, aes(x = Description, y = log10p)) +
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
p


