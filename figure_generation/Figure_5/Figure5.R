#==============================================#
# LRS WGS QC #
# Figure-5#
#==============================================#
library(data.table)
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggpubr)

setwd("/media/london_A/mengxin/GTOP_code/fig-5")

# Fig.5a ------------------------------------------------------------------
# load data
df_sv_eqtl.r2 <- readRDS("./input/dat_sv.RDS")
df_tr_eqtl.r2 <- readRDS("./input/dat_tr.RDS")
df_tissue_color <- readRDS("./input/tissue_colors.RDS")

# poorly tagged eSV: LD R2<0.8
df_sv.poor <- df_sv_eqtl.r2 %>% filter(R2<0.8);nrow(df_sv.poor)/nrow(df_sv_eqtl.r2) # 23.7%
df_tr.poor <- df_tr_eqtl.r2 %>% filter(R2<0.8);nrow(df_tr.poor)/nrow(df_tr_eqtl.r2) # 33.4%

# no. of poorly tagged SV-eQTL and TR-eQTL
median(table(df_sv.poor$Tissue))
median(table(df_tr.poor$Tissue))

df_sv_eqtl.r2$R2_pseu <- df_sv_eqtl.r2$R2
df_sv_eqtl.r2$R2_pseu[df_sv_eqtl.r2$R2<0.2] <- 0.2
df_tr_eqtl.r2$R2_pseu <- df_tr_eqtl.r2$R2
df_tr_eqtl.r2$R2_pseu[df_tr_eqtl.r2$R2<0.2] <- 0.2


p1 <- ggplot(df_sv_eqtl.r2,aes(x=R2_pseu,fill=Tissue)) + geom_histogram(position = "stack") + theme_pubr() + 
  scale_fill_manual(breaks = df_tissue_color$Tissue,values = paste0("#",df_tissue_color$Tissue_Color_Code)) + 
  scale_x_continuous(breaks = c(0.2,0.4,0.6,0.8,1.0));p1
p2 <- ggplot(df_tr_eqtl.r2,aes(x=R2_pseu,fill=Tissue)) + geom_histogram(position = "stack") + theme_pubr() + 
  scale_fill_manual(breaks = df_tissue_color$Tissue,values = paste0("#",df_tissue_color$Tissue_Color_Code)) + 
  theme(legend.position = "none") + scale_x_continuous(breaks = c(0.2,0.4,0.6,0.8,1.0));p2
p<-cowplot::plot_grid(p1+theme(legend.position = "none"),p2,ncol=2,align = "v")

p

# Fig.5b ------------------------------------------------------------------


df_plot_count <- readRDS("./input/joint_finemap.res.RDS")


p3 <- ggplot(df_plot_count,aes(x=Var1,y=log2(Freq))) + geom_bar(stat = "identity",width=.6,aes(fill = Var2),position = position_dodge()) + theme_pubr() + 
  scale_fill_manual(breaks = c("SNV","TR","SV"),values = c("#8090B4","#963628","#227E84")) + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1, vjust = .5));p3


# Fig.5c ------------------------------------------------------------------


df_plot <- readRDS("./input/fig.c1.input.RDS")
df_plot.wl <- readRDS("./input/fig.c2.input.RDS")

library(ggplot2)
library(ggpubr)


p1 <- ggplot(df_plot,aes(x=Var1,y=Freq,fill=Var2)) + geom_bar(stat = "identity",position = "stack") + theme_pubr() + 
  scale_fill_manual(breaks = c(0,1,2),values = c("#8090B4","#8090B4","#8090B4")) + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5)) + 
  xlab("Tissues") + ylab("# of credible sets") + theme(legend.position = "none");p1

library(reshape2)

p2 <- ggplot(df_plot.wl,aes(x=Tissue,y=value,fill=variable)) + geom_bar(stat = "identity",position = "stack") + theme_pubr() + 
  scale_fill_manual(breaks = c("PC_snv","PC_sv_tr","PC_sv_trLead"),values = c("#bdbdbd","#fc9272","#de2d26")) + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5)) + 
  xlab("Tissues") + ylab("# of credible sets");p2

cowplot::plot_grid(p1,p2,ncol = 1,align = "v")


# Fig.5d ------------------------------------------------------------------


df_comb <- readRDS("./input/Fig5d_data.RDS")

df_comb$cs_class <- factor(df_comb$cs_class,levels = c("SNV_TR","SNV_only","TR_only","SNV_SV_TR","SNV_SV","SV_only","SV_TR"))
df_comb.snv <- df_comb %>% filter(cs_class =="SNV_only")
df_comb.sv <- df_comb %>% filter(cs_class %in% c("SNV_SV_TR","SNV_SV","SV_only","SV_TR"))
df_comb.tr <- df_comb %>% filter(cs_class %in% c("SNV_SV_TR","SNV_TR","TR_only","SV_TR"))
df_comb.sv %<>% group_by(Var1) %>% mutate(Count=sum(Freq)) %>% select(Var1,Count) %>% distinct()
df_comb.sv$cs_class <- "cs_sv"
df_comb.tr %<>% group_by(Var1) %>% mutate(Count=sum(Freq)) %>% select(Var1,Count) %>% distinct()
df_comb.tr$cs_class <- "cs_tr"
names(df_comb.snv) <- c("Var1","Count","cs_class")
df_plot <- rbind(df_comb.snv,df_comb.sv,df_comb.tr)
df_plot$Var1 <- as.character(df_plot$Var1)
df_plot$cs_class <- as.character(df_plot$cs_class)
df_plot$Var1[df_plot$Var1=="SNV" & df_plot$cs_class=="SNV_only"] <- "not_lead"
df_plot$Var1[df_plot$Var1=="SV" & df_plot$cs_class=="cs_sv"] <- "lead"
df_plot$Var1[df_plot$Var1=="TR" & df_plot$cs_class=="cs_tr"] <- "lead"

df_plot$Var1[df_plot$Var1!="lead" & df_plot$cs_class=="cs_sv"] <- "not_lead"
df_plot$Var1[df_plot$Var1!="lead" & df_plot$cs_class=="cs_tr"] <- "not_lead"
df_plot$Var1 <- factor(df_plot$Var1,levels = c("not_lead","lead"))


p6 <- ggplot(df_plot,aes(x=cs_class,y=Count,fill=Var1)) + geom_bar(stat = "identity",position = "stack",width=.85) + theme_pubr() + 
  scale_fill_manual(breaks = c("not_lead","lead"),values = c("#bdbdbd","#FC9272")) +
  xlab("CS groups") + ylab("# of credible sets");p6


# Fig.5e ------------------------------------------------------------------

# load tissue color
df_tissue_color <- readRDS("./input/tissue_colors.RDS")

df_fm <- readRDS("./input/SV_enrich_causal_CS.compare_to_SNV.RDS")
df_fm$PIP_num <- df_fm$PIP
df_fm$PIP <- factor(df_fm$PIP)

df_tissue_color %<>% filter(Tissue %in% unique(df_fm$Tissue))

p1 <- ggplot(df_fm,aes(x=PIP,y=OR)) + geom_boxplot(outlier.shape = NA,aes(fill = PIP_num),alpha=.8) +
  geom_jitter(aes(color=Tissue),width=.2,size=1) + theme_pubr() +
  scale_color_manual(breaks = df_tissue_color$Tissue,values = paste0("#",df_tissue_color$Tissue_Color_Code)) +
  scale_fill_gradient(low = "#fff7ec",high="#b30000")
print(p1)

res <- lm(OR ~ PIP_num, data=df_fm)
summary(res)

df_m <- readRDS("./input/TR_enrich_causal_CS.compare_to_SNV.RDS")
df_m$PIP_num <- df_m$PIP
df_m$PIP <- factor(df_m$PIP)
df_tissue_color %<>% filter(Tissue %in% unique(df_m$Tissue))

p2 <- ggplot(df_m,aes(x=PIP,y=OR)) + geom_boxplot(outlier.shape = NA,aes(fill = PIP_num),alpha=.8) +
  geom_jitter(aes(color=Tissue),width=.2,size=1) + theme_pubr() +
  scale_color_manual(breaks = df_tissue_color$Tissue,values = paste0("#",df_tissue_color$Tissue_Color_Code)) +
  scale_fill_gradient(low = "#fff7ec",high="#b30000");p2
print(p2)

res <- lm(OR ~ PIP_num, data=df_m)
summary(res)
cowplot::plot_grid(p1,p2,ncol=2)
