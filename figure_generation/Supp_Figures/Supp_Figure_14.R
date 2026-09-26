#===============================================#
# MS peptide  #
# Supp-Figure-14 #
#===============================================#
library(ggplot2)
library(ggpubr)
library(tidyverse)

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig14. Figure S10/input")

# Supp.Fig.14a -peptide length distribution--------------------------------

theme_gtop <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_blank()
    )
}
data <- read.delim( "Fig_S14a.txt", check.names = FALSE)
ggplot(data, aes(peptide_length, count)) +
  geom_col(fill = "#4C72B0", width = 0.85) +
  labs(x = "Peptide length (aa)", y = "Number of peptides") +
  theme_gtop()+
  scale_y_continuous(breaks = seq(0,13500,by=1500))

# Supp.Fig.14b tissue-level peptide support for annotated transcripts----------

library(dplyr)
theme_gtop <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_blank()
    )
}
pretty_label <- function(x) {
  x <- gsub("_", " ", as.character(x))
  x <- tolower(x)
  ifelse(nchar(x) > 0, paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x))), x)
}
first_existing_col <- function(data, cols, default = NA_character_) {
  for (col in cols) {
    if (col %in% names(data)) return(data[[col]])
  }
  rep(default, nrow(data))
}
data <- read.delim("Fig_S14b.txt", check.names = FALSE)
data <- data[data$validate_status %in% c("Unique peptide", "Shared peptide"), ]
data$tissue_label <- first_existing_col(data, c("tissue_label"), NA_character_)
data$tissue_label[is.na(data$tissue_label)] <- pretty_label(data$tissue[is.na(data$tissue_label)])
data$tissue_color <- first_existing_col(data, c("tissue_color"), "#7f8080")
data$tissue_label <- factor(data$tissue_label, levels = sort(unique(as.character(data$tissue_label))))
data$validate_status <- factor(data$validate_status, levels = c("Unique peptide", "Shared peptide"))
tissue_points <- data[!duplicated(data$tissue_label), c("tissue_label", "tissue_color")]
ggplot(data, aes(tissue_label, ratio, fill = validate_status)) +
  geom_col(width = 0.8, position = position_stack(reverse = TRUE)) +
  geom_point(
    data = tissue_points,
    aes(x = tissue_label, y = -0.03, color = tissue_color),
    inherit.aes = FALSE,
    size = 2
  ) +
  scale_fill_manual(values = c("Unique peptide" = "#5b8fc7", "Shared peptide" = "#a3c1e1")) +
  scale_color_identity() +
  scale_y_continuous(breaks = seq(0,1,by=0.2))+
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  labs(x = "Tissues", y = "Proportion of isoform with peptides") +
  theme_gtop() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Supp.Fig.14c -unique-peptide support count groups for novel isoforms--------------------

theme_gtop <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_blank()
    )
}
pretty_label <- function(x) {
  x <- gsub("_", " ", as.character(x))
  x <- tolower(x)
  ifelse(nchar(x) > 0, paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x))), x)
}
data <- read.delim("Fig_S14c.txt", check.names = FALSE)
data$peptide_group <- factor(data$peptide_group, levels = c("1", "2 - 5", "> 5"))
data$tissue <- pretty_label(data$tissue)
ggplot(data, aes(tissue, count, fill = peptide_group)) +
  geom_col(width = 0.8, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = c("1" = "#deebf7", "2 - 5" = "#9ecae1", "> 5" = "#3182bd"), drop = FALSE) +
  labs(x = NULL, y = "Number of transcripts") +
  theme_gtop()


# Supp.Fig.14d ------------------------------------------------------------
tsv <- "02_structure_pairs_with_categories.tsv"

d   <- read.delim(tsv, stringsAsFactors = FALSE)
chg <- d$structure_status == "Structure changed"          # TM-score < 0.5

sel <- list(all  = rep(TRUE, nrow(d)),
            gain = d$domain_category == "Domain gain",
            loss = d$domain_category == "Domain loss")
pct <- sapply(sel, function(i) 100 * sum(chg[i]) / sum(i))  # 51.4 / 43.5 / 85.6

COL <- c(no_change = "#9E9E9F", changed = "#246ca0")
x   <- c(all = 0.80, gain = 3.20, loss = 4.33)              # bar centres
w   <- 0.78                                                 # bar width

panel_d <- function() {
  par(mai = c(0.90, 1.29, 0.32, 0.10), xpd = NA, family = "Arial")
  plot(NA, xlim = c(0, 9.43), ylim = c(0, 100), axes = FALSE, xlab = "", ylab = "",
       xaxs = "i", yaxs = "i")
  axis(2, at = seq(0, 100, 25), labels = paste0(seq(0, 100, 25), "%"),
       las = 1, cex.axis = 1.0)                             # y axis (left)
  mtext("Percentage", side = 2, line = 2.6, cex = 1.05)
  segments(c(0, 2.5), 0, c(1.6, 5.1), 0)                    # x axis, split in two groups
  segments(2.5, 0, 2.5, 100)                                # divider between the groups
  
  rect(x - w/2, 0,   x + w/2, pct, col = COL["changed"],   border = "black")
  rect(x - w/2, pct, x + w/2, 100, col = COL["no_change"], border = "black")
  
  text(x, -7, c("all", "gain", "loss"), font = 2, cex = 1.05)
  text(0.80,  -19.5, "Tested transcripts",    cex = 1.05)
  text(3.765, -19.5, "Protein domain change", cex = 1.05)
  text(-0.35, 106, "d", font = 2, cex = 1.5)
  
  legend(5.13, 102, xjust = 0, yjust = 1, bty = "n", border = "black",
         cex = 1.0, y.intersp = 1.25,
         fill = COL[c("no_change", "changed")],
         legend = c(sprintf("No structure change:\nTM-score \u2265 0.5 (n=%d)", sum(!chg)),
                    sprintf("Structure changed:\nTM-score < 0.5 (n=%d)",     sum(chg))))
}

panel_d()
