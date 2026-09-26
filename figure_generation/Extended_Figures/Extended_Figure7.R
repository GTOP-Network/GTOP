#==================================#
# Extended Figure-7
#==================================#

setwd("/media/london_A/mengxin/GTOP_code/extend/extend_7")

library(data.table)
library(tidyverse)
library(reticulate)


## Extended.Fig.7a, fd-sQTLs --------------------------
library(reticulate)
py_require("matplotlib==3.6.3")
py_require("git+https://github.com/aabiddanda/geovar")
py_config()

py_run_string(
  '
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# New NumPy removed np.row_stack
if not hasattr(np, "row_stack"):
    np.row_stack = np.vstack

# New Matplotlib removed matplotlib.cm.get_cmap
if not hasattr(cm, "get_cmap"):
    def get_cmap(name=None, lut=None):
        cmap = plt.colormaps[name if name is not None else "viridis"]
        if lut is not None:
            cmap = cmap.resampled(lut)
        return cmap
    cm.get_cmap = get_cmap

# Import geovar after compatibility patches
from geovar import *

plt.rcParams["pdf.fonttype"] = 42

geovar_test = GeoVar()
geovar_test.add_freq_mat("input/ext_Fig7a.txt")
geovar_test.geovar_binning()

geovar_plot = GeoVarPlot()
geovar_plot.add_data_geovar(geovar_test)
geovar_plot.filter_data()
geovar_plot.add_cmap()

fig, ax = plt.subplots(1, 1, figsize=(3, 6))
geovar_plot.plot_geovar(ax)

ax.set_xticklabels(geovar_plot.poplist)

plt.savefig(
    "fd_sQTL_freq_data.pdf",
    dpi=600,
    bbox_inches="tight"
)

plt.close()
'
)


## Extended Data Fig.7b mash of eQTL portability ---------------------------------------------
color_vec <- readRDS("../../fig-4/input/tissue_color.RDS")

count_df <- fread("./input/ext_Fig7b.txt")
count_df$type2 <- factor(count_df$type2, levels = c("nominal", "mash"))

ggplot(
  count_df,
  aes(x = type2, y = ratio, fill = type2)
) +
  geom_boxplot(outlier.color = "NA") +
  geom_line(
    aes(group = tissue, color = tissue),
    position = position_dodge(width = 0.4),
    alpha = 0.5
  ) +
  geom_point(
    aes(group = tissue, color = tissue),
    position = position_dodge(width = 0.4)
  ) +
  # ggbeeswarm::geom_quasirandom(aes(color = tissue), size = 2) +
  scale_color_manual(values = color_vec) +
  theme_classic() +
  labs(x = "", y = "Proportion of portable eGenes") +
  ggpubr::stat_compare_means(
    method = "t.test",
    paired = T
  ) +
  scale_fill_manual(values = c("#6874b4", "#4bb9b9"))
