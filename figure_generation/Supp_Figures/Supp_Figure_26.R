#====================================#
# join finemapping res of sQTL #
# # Supp-Figure-26 # #
#===================================#

library(data.table)
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggpubr)

setwd("/path/to/GTOP_code/supp/supp_fig26")
df_plotall <- fread("./input/Figure S26.txt")

# supp.Figure.26a join finemapping res of juQTL ---------------------------


df_plot <- df_plotall %>% dplyr::filter(qtltype=="juQTL")
order <- c(
  "Adrenal_Gland","Adipose",  "Pancreas_Tail", "Pancreas_Head","Pancreas_Body",
  "Liver", "Gallbladder",  "Whole_Blood", "Spleen", "Muscle", "Skin"
)
df_plot$Var1<-factor(df_plot$Var1,levels = order )
p1 <- ggplot(df_plot,aes(x=Var1,y=Freq,fill=Var2)) + geom_bar(stat = "identity",position = "stack", fill="#8090b4") +
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5)) + 
  xlab("Tissues") + ylab("# of credible sets");p1

library(reshape2)
df_plot.w <- dcast(df_plot,Var1 ~ Var2, value.var = "Freq")
names(df_plot.w) <- c("Tissue","SNV","SV_TR","SV_TR_lead")
df_plot.w$total <- apply(df_plot.w[,-1],1,sum)

df_plot.w %<>% mutate(PC_snv = SNV/total, PC_sv_tr=SV_TR/total, PC_sv_trLead=SV_TR_lead/total) %>% select(Tissue,PC_snv,PC_sv_tr,PC_sv_trLead)

df_plot.wl <- melt(df_plot.w,id.vars = "Tissue")
df_plot.wl$Tissue <- as.character(df_plot.wl$Tissue)
df_plot.wl$Tissue <- factor(df_plot.wl$Tissue, levels = order)
p2 <- ggplot(df_plot.wl,aes(x=Tissue,y=value,fill=variable)) + geom_bar(stat = "identity",position = "stack") + theme_pubr() + 
  scale_fill_manual(breaks = c("PC_snv","PC_sv_tr","PC_sv_trLead"),values = c("#bdbdbd","#fc9272","#de2d26")) + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5)) + 
  xlab("Tissues") + ylab("# of credible sets");p2



# supp.Figure.26b join finemapping res of tuQTL ---------------------------


df_plot <- df_plotall %>% dplyr::filter(qtltype=="tuQTL")
order <- c(
  "Adrenal_Gland","Adipose","Liver","Pancreas_Head","Pancreas_Tail","Pancreas_Body",
   "Gallbladder",  "Whole_Blood", "Spleen", "Muscle", "Skin"
)
df_plot$Var1<-factor(df_plot$Var1,levels = order )
p3 <- ggplot(df_plot,aes(x=Var1,y=Freq)) + geom_bar(stat = "identity",position = "stack", fill="#8090b4") + 
  theme_pubr() + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5)) + 
  xlab("Tissues") + ylab("# of credible sets");p3

library(reshape2)
df_plot.w <- dcast(df_plot,Var1 ~ Var2, value.var = "Freq")
names(df_plot.w) <- c("Tissue","SNV","SV_TR","SV_TR_lead")
df_plot.w$total <- apply(df_plot.w[,-1],1,sum)

df_plot.w %<>% mutate(PC_snv = SNV/total, PC_sv_tr=SV_TR/total, PC_sv_trLead=SV_TR_lead/total) %>% select(Tissue,PC_snv,PC_sv_tr,PC_sv_trLead)
df_plot.wl <- melt(df_plot.w,id.vars = "Tissue")
df_plot.wl$Tissue <- as.character(df_plot.wl$Tissue)
df_plot.wl$Tissue <- factor(df_plot.wl$Tissue, levels = order)
p4 <- ggplot(df_plot.wl,aes(x=Tissue,y=value,fill=variable)) + geom_bar(stat = "identity",position = "stack") + theme_pubr() + 
  scale_fill_manual(breaks = c("PC_snv","PC_sv_tr","PC_sv_trLead"),values = c("#bdbdbd","#fc9272","#de2d26")) + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = .5)) + 
  xlab("Tissues") + ylab("# of credible sets");p4



