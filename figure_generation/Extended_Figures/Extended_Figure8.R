#==================================#
# SV-sQTLs and TR-sQTLs #
# Extended Figure-8 #
#==================================#

library(ggplot2)
library(ggpubr)

setwd("/media/london_A/mengxin/GTOP_code/extend/extend_8")


# Extended.Fig.8 ----------------------------------------------------------

df_tissue_color <- readRDS("./input/tissue_colors.RDS")

df_sv_sqtl.r2 <- readRDS("./input/dat_sv.sQTL.RDS")
df_sv_sqtl.r2$R2_pseu <- df_sv_sqtl.r2$R2
df_sv_sqtl.r2$R2_pseu[df_sv_sqtl.r2$R2<0.2] <- 0.2
p1 <- ggplot(df_sv_sqtl.r2,aes(x=R2_pseu,fill=Tissue)) + geom_histogram(position = "stack") + theme_pubr() + 
  scale_fill_manual(breaks = df_tissue_color$Tissue,values = paste0("#",df_tissue_color$Tissue_Color_Code)) + 
  theme(legend.position = "none") + scale_x_continuous(breaks = c(0.2,0.4,0.6,0.8,1.0));p1


df_tr_sqtl.r2 <- readRDS("./input/dat_tr.sQTL.RDS")
df_tr_sqtl.r2$R2_pseu <- df_tr_sqtl.r2$R2
df_tr_sqtl.r2$R2_pseu[df_tr_sqtl.r2$R2<0.2] <- 0.2
p2 <- ggplot(df_tr_sqtl.r2,aes(x=R2_pseu,fill=Tissue)) + geom_histogram(position = "stack") + theme_pubr() + 
  scale_fill_manual(breaks = df_tissue_color$Tissue,values = paste0("#",df_tissue_color$Tissue_Color_Code)) + 
  theme(legend.position = "none") + scale_x_continuous(breaks = c(0.2,0.4,0.6,0.8,1.0));p2
cowplot::plot_grid(p1+theme(legend.position = "none"),p2,ncol=1,align = "v")
