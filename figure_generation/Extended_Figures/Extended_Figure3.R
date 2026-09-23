#==============================================#
# ASE #
# Extended-Figure-3#
#==============================================#


library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)

setwd("/media/london_A/mengxin/GTOP_code/extend/extend_3")
# Extended.Fig.3ab ASE & ASTS ----------------------------------------------
df_plot <- fread("./input/Exfig 3a.txt")
reg <- lm(formula = sig ~ total,data=df_plot)
coeff <- coefficients(reg)
intercept <- coeff[1]
slope <- coeff[2]
p1 <- ggplot(df_plot,aes(x=total,y=sig,color=Abbreviation)) + geom_point(size=3) + theme_pubr() + 
  scale_color_manual(breaks = df_plot$Abbreviation,values = paste0("#",df_plot$Tissue_Color_Code)) + 
  xlab("Number of genes tested") + ylab("Number of significant genes") + 
  geom_abline(intercept = intercept,slope = slope,color="black",linetype="dashed",size=1.5)
p1

# Extended.fig.3b ---------------------------------------------------------

df_plot <- fread("./input/Exfig 3b.txt")
reg <- lm(formula = sig.asts ~ total.asts,data=df_plot)
coeff <- coefficients(reg)
intercept <- coeff[1]
slope <- coeff[2]
p2 <- ggplot(df_plot,aes(x=total.asts,y=sig.asts,color=Abbreviation)) + geom_point(size=3) + theme_pubr() + 
  scale_color_manual(breaks = df_plot$Abbreviation,values = paste0("#",df_plot$Tissue_Color_Code)) + 
  xlab("Number of genes tested") + ylab("Number of significant genes") + 
  geom_abline(intercept = intercept,slope = slope,color="black",linetype="dashed",size=1.5)

cowplot::plot_grid(p1,p2,align = "h",ncol=2)


# Extended.fig.3c  ---------------------------------------------------------
df_plot <- fread("./input/Exfig 3c.txt")
ggplot(df_plot,aes(x=Var1,y=Freq)) + 
  geom_bar(stat = "identity",width=.6, fill="#6c81b1") + 
  xlab("# transcripts per gene")+
  ylab("# Genes")+
  theme_pubr()


# Extended.fig.3d  ---------------------------------------------------------

df_sample.l <- fread("./input/Exfig 3d.txt") %>% 
  mutate(Var2=as.character(Var2))
p3 <- ggplot(df_sample.l,aes(x=Var1,y=Freq,fill=Var2)) + geom_bar(stat = "identity",width=.8) + theme_pubr() + 
  theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=.5)) + xlab("")+ylab("# ASE events")+
  scale_y_log10() + 
  scale_fill_manual(breaks = c("0","1"),labels=c("non-significant events", "significant events"),values = c("grey","#4380b8"))+
  theme(axis.text.x = element_blank());p3

# Extended.fig.3e  ---------------------------------------------------------
df_sample <- fread("./input/Exfig 3e.txt") %>% 
  mutate(Var2=as.character(Var2))
p4 <- ggplot(df_sample,aes(x=Var1,y=Freq,fill=Var2)) + geom_bar(stat = "identity",width=.8) + theme_pubr() + 
  theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=.5)) + xlab("")+ylab("# ASTS events")+
  scale_y_log10() + 
  scale_fill_manual(breaks = c("0","1"), labels=c("non-significant events", "significant events"),values = c("grey","#4380b8"))+
  theme(axis.text.x = element_blank());p4






