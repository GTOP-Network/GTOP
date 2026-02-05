#==============================================#
# Short-reads-RNA-QC #
# Supp-Figure-12#
#==============================================#
setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig3/input")

library(data.table)
library(readxl)
library(ggplot2)
qc_data <- read_excel("RNA_matrix_1586.xlsx")

colnames(qc_data)

# Supp.Figure.12d Distribution of Exonic Rate -----------------------------


library(scales) 

ggplot(qc_data, aes(x = `Exonic Rate` * 100)) +
  geom_density(aes(y = ..density..), fill = "grey80") +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),expand = c(0, 0)) +  
  scale_y_continuous(expand = c(0, 0)) +  
  labs(x = "Exonic Rate (%)", y = "Density",
       title="Distribution of Exonic Rate") +
  theme_classic()+
  theme(
    plot.title = element_text(hjust = 0.5)  
  )


# Supp.Figure.12e Distribution of Intronic Rate ---------------------------


ggplot(qc_data, aes(x = `Intronic Rate` * 100)) +
  geom_density(aes(y = ..density..), fill = "grey80") +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),expand = c(0, 0)) +  
  scale_y_continuous(expand = c(0, 0)) +  
  labs(x = "Intronic Rate (%)", y = "Density",
       title="Distribution of Intronic Rate") +
  theme_classic()+
  theme(
    plot.title = element_text(hjust = 0.5)  
  )


# Supp.Figure.12f Distribution of Intergenic Rate --------------------------

ggplot(qc_data, aes(x = `Intergenic Rate` * 100)) +
  geom_density(aes(y = ..density..), fill = "grey80") +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),expand = c(0, 0)) +  
  scale_y_continuous(expand = c(0, 0)) +  
  labs(x = "Intergenic Rate (%)", y = "Density",
       title="Distribution of Intergenic Rate") +
  theme_classic()+
  theme(
    plot.title = element_text(hjust = 0.5)  
  )


# Supp.Figure.12g Distribution of Intragenic Rate -------------------------

ggplot(qc_data, aes(x = `Intragenic Rate` * 100)) +
  geom_density(aes(y = ..density..), fill = "grey80") +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),expand = c(0, 0)) +  
  scale_y_continuous(limits = c(0, 0.15), breaks = seq(0, 0.15, 0.05),expand = c(0, 0)) + 
  labs(x = "Intragenic Rate (%)", y = "Density",
       title="Distribution of Intragenic Rate") +
  theme_classic()+
  theme(
    plot.title = element_text(hjust = 0.5)  
  )


# Supp.Figure.12h Distribution of rRNA Rate (log10-transformed) ------------

ggplot(qc_data, aes(x = `rRNA Rate`*100)) +
  geom_density(fill = "grey80") +
  scale_x_continuous(
    trans = "log10",
    breaks = c(1e-4, 1e-2, 1),
    labels = c("0.0001", "0.01", "1"),
    expand = c(0, 0)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x = "rRNA Rate (%)",
    y = "Density",
    title = "Distribution of rRNA Rate"
  ) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))



# Supp.Figure.12i Distribution of Median 3' Bias ---------------------------

ggplot(qc_data, aes(x = `Median 3' bias`)) +
  geom_density(aes(y = ..density..), fill = "grey80") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),expand = c(0, 0)) +  
  scale_y_continuous(expand = c(0, 0)) + 
  labs(x = "Median 3' bias", y = "Density",
       title="Distribution of Median 3' Bias") +
  theme_classic()+
  theme(
    plot.title = element_text(hjust = 0.5)  
  )


# Supp.Figure.12j Mapped reads-------------------------------------------


library(ggplot2)
library(dplyr)
library(stringr)
library(data.table)
load("map_reads_plot.RData")

#set order
load("order.RData")
length(unique(plot_data$Tissue))
order<-unique(order[order$order %in% unique(plot_data$Tissue),,drop=F])

plot_data$Tissue <- factor(plot_data$Tissue, levels = order$order)

tissue_colors <- plot_data %>%
  distinct(Tissue, Tissue_Color_Code) %>%
  pull(Tissue_Color_Code)

names(tissue_colors) <- unique(plot_data$Tissue)



# violin -----------------------------------------------------------------

p <- ggplot(plot_data, aes(x = Tissue, y = total_reads, fill = Tissue, alpha = Dataset)) +
  geom_violin(trim = FALSE, scale = "width", 
              color = "black", linewidth = 0.3, adjust = 1.2,
              position = position_dodge(width = 0.8)) +
  geom_boxplot(width = 0.12, alpha = 0.8, 
               color = "black", linewidth = 0.2,
               position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 82, color = "blue", linetype = "dashed", linewidth = 0.8, alpha = 0.6) +
  geom_hline(yintercept = 193.365406, color = "blue", linetype = "dashed", linewidth = 0.8, alpha = 0.7) +
  
  scale_fill_manual(values = tissue_colors) +
  scale_alpha_manual(values = c("gtex" = 0.4, "GTOP" = 1)) +  
  labs(y = "Mapped reads (millions)", fill = "Tissue") +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    panel.grid = element_blank()
  )

print(p)




