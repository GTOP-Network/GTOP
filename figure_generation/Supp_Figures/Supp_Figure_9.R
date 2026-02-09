#===============================================#
# LRS read QC #
# Supp-Figure-9                         #
#===============================================#
library(ggplot2)
library(ggpubr)
library(ComplexUpset)
library(tidyverse)
library(data.table)

setwd("/media/london_A/mengxin/GTOP_code/supp/supp_fig9/input")


# Supp.Fig.9a Candidate transcripts --------------------------------------------

df <- read_tsv("supp9a.LR_candidate_isoform_tool_upset.txt",
  show_col_types = FALSE
)
df <- df %>%
  rename(
    Bambu   = bambu,
    FLAIR   = flair,
    `Iso-Seq` = isoseq
  )
tool_cols <- c("Bambu", "FLAIR", "Iso-Seq")

upset(
  df,
  intersect = tool_cols,
  base_annotations = list(
    "Intersection size" = intersection_size(
      counts  = FALSE
    )
  ),
  
  width_ratio = 0.3
)  +
  theme(
    text = element_text(size = 12, color = "black"),
    axis.text  = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.line  = element_line(color = "black", linewidth = 0.6),
    panel.background = element_blank(),
    panel.grid = element_blank(), 
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11),
    strip.background = element_blank(),
    strip.text = element_text(color = "black")
  )


# Supp.Fig.9c FLNC read counts --------------------------------------------

# Load data
df <- fread("supp9c.LR_isoform_QC_reads.txt.gz")
flnc_reads <- df[, rowSums(.SD), .SDcols = patterns("^")]
log_flnc_reads <- log10(flnc_reads)

original_threshold <- 10
log_threshold <- log10(original_threshold)  # exactly 1.0

# Create bins
max_log_val <- ceiling(max(log_flnc_reads, na.rm = TRUE) * 10) / 5
custom_bins <- seq(0, max_log_val + 0.2, by = 0.2)

bin_table <- tibble(
  bin_left  = custom_bins[-length(custom_bins)],
  bin_right = custom_bins[-1]
)

binned_counts <- tibble(log_flnc = log_flnc_reads) %>%
  mutate(
    bin = cut(
      log_flnc,
      breaks = custom_bins,
      include.lowest = TRUE,
      right = FALSE
    )
  ) %>%
  count(bin)

plot_df <- bin_table %>%
  mutate(
    bin = factor(
      paste0("[", bin_left, ",", bin_right, ")"),
      levels = levels(binned_counts$bin)
    )
  ) %>%
  left_join(binned_counts, by = "bin") %>%
  mutate(
    n = ifelse(is.na(n), 0, n),
    pass = ifelse(bin_left >= log_threshold, "passed", "filtered")
  )
plot_df$pass <- factor(plot_df$pass, levels = c("filtered", "passed"))
legend_labels <- c(
  "filtered" = paste0("Filtered out (< ", original_threshold, ")"),
  "passed"   = paste0("Passed filters (≥ ", original_threshold, ")")
)

# Plot
ggplot(plot_df) +
  geom_rect(
    aes(
      xmin = bin_left,
      xmax = bin_right,
      ymin = 0,
      ymax = n,
      fill = pass
    ),
    color = "black",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = log_threshold,
    linetype = "dashed",
    linewidth = 1,
    color = "red"
  ) +
  scale_fill_manual(
    values = c(
      "filtered" = "#FFA07A",
      "passed"   = "#4682B4"
    ),
    labels = legend_labels,
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(0,log10(max(flnc_reads))+0.2),
    breaks = seq(0, max_log_val + 1, by = 1),
    expand = c(0, 0)
  ) +
  labs(
    x = "Supported FLNC read count (log10 scale)",
    y = "Number of transcripts",
    fill = NULL
  ) +
  theme_pubr() +
  theme(
    legend.position = "top"
  )


# Supp.Fig.9d -------------------------------------------------------------


plot_df <- fread("supp9d.LR_isoform_QC_support_samples.txt.gz")
binwidth <- (max(plot_df$support_samples) - min(plot_df$support_samples)) / 20
ggplot(plot_df, aes(x = support_samples)) +
  geom_histogram(binwidth = binwidth, 
                 fill = "#7f8080", center = min(plot_df$support_samples) + binwidth/2,
                 color = "black") +
  xlab("Number of supporting samples") +
  ylab("Frequency") +
  theme_pubr()



# Supp.Fig.9e -------------------------------------------------------------


plot_df <- fread("supp9e.LR_read_QC_intra_priming.txt")


ggplot(plot_df) +
  geom_bar(
    aes(x = A_number, y = frequency, fill = type),
    stat = "identity",
    color = "black",
    size = 0.2, 
    width = 1     
  ) +
  scale_fill_manual(values = setNames(plot_df$color, plot_df$type)) +
  xlab("# A in 20bp downstream TTS") +
  ylab("Number of transcripts") +
  theme_pubr() +
  theme(legend.position = "none")


# Supp.Fig.9f -------------------------------------------------------------


plot_df <- fread("supp9f.LR_read_QC_pass_category.Passed filters.txt")
color <- plot_df %>% distinct(type, color)

plot_df$category <- factor(plot_df$category, levels = plot_df$category)

ggplot(plot_df)+
  geom_bar(aes(x=category, y=count, fill=type), stat = "identity", color = "black")+
  scale_fill_manual(values=setNames(color$color, color$type))+
  xlab("")+
  ylab("Transcript count")+
  theme_pubr()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


# Supp.Fig.9g -------------------------------------------------------------


plot_df <- fread("supp9g.LR_read_QC_pass_category.Filtered out.txt")
color <- plot_df %>% distinct(type, color)
plot_df$category <- factor(plot_df$category, levels = plot_df$category)

ggplot(plot_df)+
  geom_bar(aes(x=category, y=count, fill=type), stat = "identity", color = "black")+
  scale_fill_manual(values=setNames(color$color, color$type))+
  xlab("")+
  ylab("Transcript count")+
  theme_pubr()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


# Supp.Fig.9h -------------------------------------------------------------

df <- read_tsv(
  "supp9h.LR_final_isoform_tool_with_cate_upset.txt",
  show_col_types = FALSE
)

df <- df %>%
  rename(
    Bambu   = bambu,
    FLAIR   = flair,
    `Iso-Seq` = isoseq
  )
tool_cols <- c("Bambu", "FLAIR", "Iso-Seq")

sorted_cate <- c('FSM', 'ISM', 'NIC', 'NNC', 'Fusion', 'Antisense', 'Genic', 'Intergenic')

df <- df %>%
  mutate(
    struct_class = factor(
      struct_class,
      levels = c(
        intersect(sorted_cate, unique(struct_class)),
        setdiff(unique(struct_class), sorted_cate)
      )
    )
  )

fill_colors <- c(
  "FSM"       = "#226ca0",
  "ISM"       = "#f39800",
  "NIC"       = "#65a066",
  "NNC"       = "#00a29a",
  "Fusion"    = "#595757",
  "Antisense" = "#e94651",
  "Genic"     = "#7861a8",
  "Intergenic"= "#c17272"
)


upset(
  df,tool_cols,
  base_annotations=list(
    'Intersection size'=intersection_size(
      counts=FALSE,
      mapping=aes(fill=struct_class)
    ) + scale_fill_manual(values=fill_colors)
  ),
  width_ratio=0.3
)  +
  theme(
    text = element_text(size = 12, color = "black"),
    axis.text  = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.line  = element_line(color = "black", linewidth = 0.6),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11),
    strip.background = element_blank(),
    strip.text = element_text(color = "black")
  )
