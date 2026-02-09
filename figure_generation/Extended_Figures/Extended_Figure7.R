#==============================================#
# Cross-ancestry fine-mapping #
# Extend-Data-Figure-7#
#==============================================#

library(data.table)
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggpubr)
setwd("/media/london_A/mengxin/GTOP_code/extend/extended_7")

# Extended.Data.Fig.7a ----------------------------------------------------

xy_gtex <- readRDS("./input/gtex.fm.RDS")
xy_gtop <- readRDS("./input/gtop.fm.RDS")

df_gtex <- as.data.frame(table(xy_gtex$Tissue,xy_gtex$share))
df_gtop <- as.data.frame(table(xy_gtop$Tissue,xy_gtop$share))
df_tissue <- as.data.frame(table(xy_gtop$Tissue))
df_tissue <- df_tissue[order(-df_tissue$Freq),]

df_gtex$Var1 <- factor(df_gtex$Var1,levels = df_tissue$Var1)
df_gtop$Var1 <- factor(df_gtop$Var1,levels = df_tissue$Var1)
df_gtex$Var2 <- factor(df_gtex$Var2,levels = c("specific","shared"))
df_gtop$Var2 <- factor(df_gtop$Var2,levels = c("specific","shared"))
p1 <- ggplot(df_gtex,aes(x=Var1,y=Freq,fill=Var2)) +
  geom_bar(stat="identity",position = position_stack()) + theme_pubr() +
  ylim(0,150000);p1
p2 <- ggplot(df_gtop,aes(x=Var1,y=Freq,fill=Var2)) +
  geom_bar(stat="identity",position = position_stack()) + theme_pubr() +
  ylim(0,150000);p2

cowplot::plot_grid(p2,p1,ncol=1,align="v")

df_gtex.w <- dcast(df_gtex,Var1~Var2,value.var = "Freq")
df_gtop.w <- dcast(df_gtop,Var1~Var2,value.var = "Freq")

df_gtex.w %<>% mutate(prp_s=specific/(specific+shared))
df_gtop.w %<>% mutate(prp_s=specific/(specific+shared))


df_gtex <- as.data.frame(table(xy_gtex$Tissue,xy_gtex$share))
df_gtop <- as.data.frame(table(xy_gtop$Tissue,xy_gtop$share))
df_tissue <- as.data.frame(table(xy_gtop$Tissue))
df_tissue <- df_tissue[order(-df_tissue$Freq),]

df_gtex$Var1 <- factor(df_gtex$Var1,levels = df_tissue$Var1)
df_gtop$Var1 <- factor(df_gtop$Var1,levels = df_tissue$Var1)
df_gtex$Var2 <- factor(df_gtex$Var2,levels = c("specific","shared"))
df_gtop$Var2 <- factor(df_gtop$Var2,levels = c("specific","shared"))
p1 <- ggplot(df_gtex,aes(x=Var1,y=Freq,fill=Var2)) +
  geom_bar(stat="identity",position = position_stack()) + theme_pubr() +
  ylim(0,150000);p1
p2 <- ggplot(df_gtop,aes(x=Var1,y=Freq,fill=Var2)) +
  geom_bar(stat="identity",position = position_stack()) + theme_pubr() +
  ylim(0,150000);p2

cowplot::plot_grid(p2,p1,ncol=1,align="v")

# Extended.Data.Fig.7b ----------------------------------------------------


rm_version_num <- function(x){
  return(strsplit(x,split=".",fixed = T)[[1]][1])
}

dat.gtex <- readRDS("./input/dat.fm_gtex.RDS")
dat.gtop <- readRDS("./input/dat.fm_gtop.RDS")


dat.gtex$locus_id <- sapply(dat.gtex$locus_id,rm_version_num)
dat.gtop$locus_id <- sapply(dat.gtop$locus_id,rm_version_num)

dat.gtex$Tissue[dat.gtex$Tissue=="Adipose_Visceral_Omentum"] <- "Adipose"
dat.gtex$Tissue[dat.gtex$Tissue=="Muscle_Skeletal"] <- "Muscle"
dat.gtex$Tissue[dat.gtex$Tissue=="Skin_Not_Sun_Exposed_Suprapubic"] <- "Skin"
dat.gtex$Tissue[dat.gtex$Tissue=="Pancreas"] <- "Pancreas_Body"

selected_tissues <- c("Adipose","Adrenal_Gland","Liver","Muscle","Pancreas_Body","Spleen","Whole_Blood")



dat.gtex %<>% filter(Tissue %in% selected_tissues)
dat.gtop %<>% filter(Tissue %in% selected_tissues)

dat.gtex$tissue_gene <- paste(dat.gtex$locus_id,dat.gtex$Tissue,sep=":")
dat.gtop$tissue_gene <- paste(dat.gtop$locus_id,dat.gtop$Tissue,sep=":")


# find single cs genes
df.gtex <- dat.gtex %>% select(locus_id,cs,Tissue) %>% distinct(.keep_all = T)
df.gtop <- dat.gtop %>% select(locus_id,cs,Tissue) %>% distinct(.keep_all = T)


df_count.gtex <- df.gtex %>% group_by(locus_id,Tissue) %>% summarise(cs_count=n())
df_count.gtop <- df.gtop %>% group_by(locus_id,Tissue) %>% summarise(cs_count=n())


# extract genes with single cs
df_count.gtex$tissue_gene <- paste(df_count.gtex$locus_id,df_count.gtex$Tissue,sep = ":")
df_count.gtop$tissue_gene <- paste(df_count.gtop$locus_id,df_count.gtop$Tissue,sep = ":")

x.gtex <- as.character(df_count.gtex$tissue_gene[df_count.gtex$cs_count==1])
x.gtop <- as.character(df_count.gtop$tissue_gene[df_count.gtop$cs_count==1])

overlap_tissuegenes <- intersect(x.gtex,x.gtop) # 5945

df.gtex <- dat.gtex %>% filter(tissue_gene %in% overlap_tissuegenes) %>% select(locus_id,cs,pip,cs_size,Tissue,tissue_gene) %>% 
  group_by(locus_id,cs,Tissue) %>% mutate(max_pip = max(pip)) %>% ungroup() %>% filter(pip==max_pip) %>% 
  select(locus_id,cs,max_pip,cs_size,Tissue,tissue_gene) %>% distinct(.keep_all = T)
df.gtop <- dat.gtop %>% filter(tissue_gene %in% overlap_tissuegenes) %>% select(locus_id,cs,pip,cs_size,Tissue,tissue_gene) %>% 
  group_by(locus_id,cs,Tissue) %>% mutate(max_pip = max(pip)) %>% ungroup() %>% filter(pip==max_pip) %>% 
  select(locus_id,cs,max_pip,cs_size,Tissue,tissue_gene) %>% distinct(.keep_all = T)


df_plot.gtex <- df.gtex[,c(3,6)]
df_plot.gtop <- df.gtop[,c(3,6)]
df_plot <- merge(df_plot.gtop,df_plot.gtex,by="tissue_gene")
names(df_plot) <- c("tissue_gene","maxPIP_gtop","maxPIP_gtex")

df_plot$maxPIP_all <- apply(df_plot[,2:3],1,max)



df_plot$PIP_bin_gtex <- cut(df_plot$maxPIP_gtex,breaks = c(0,0.2,0.4,0.6,0.8,1.0),labels = c("0-0.2","0.2-0.4","0.4-0.6","0.6-0.8","0.8-1.0"))
df_plot$PIP_bin_gtop <- cut(df_plot$maxPIP_gtop,breaks = c(0,0.2,0.4,0.6,0.8,1.0),labels = c("0-0.2","0.2-0.4","0.4-0.6","0.6-0.8","0.8-1.0"))

df_sank <- df_plot %>% group_by(PIP_bin_gtex,PIP_bin_gtop) %>% summarise(count=n())
df_sank$PIP_bin_gtex <- paste(df_sank$PIP_bin_gtex,"gtex",sep="-")
df_sank$PIP_bin_gtop <- paste(df_sank$PIP_bin_gtop,"gtop",sep="-")
df_sank$log2count <- log2(df_sank$count)


library(networkD3)
nodes <- data.frame(name=unique(c(df_sank$PIP_bin_gtex,df_sank$PIP_bin_gtop)))

df_sank$IDsource <- match(df_sank$PIP_bin_gtex,nodes$name)-1
df_sank$IDtarget <- match(df_sank$PIP_bin_gtop,nodes$name)-1
p2 <- sankeyNetwork(Links = df_sank, Nodes = nodes,
                    Source = "IDsource", Target = "IDtarget",
                    Value = "count", NodeID = "name", 
                    sinksRight=FALSE,fontSize = 10);p2

# Extended.Data.Fig.7c: CS size in single population fine-mapping  --------

dat.gtex <- readRDS("./input/dat.gtex.RDS")
dat.gtop <- readRDS("./input/dat.gtop.RDS")
dat.cran <- readRDS("./input/dat.cran.RDS")
# examine CS number per gene
df.gtex <- dat.gtex %>% select(locus_id,cs,Tissue) %>% distinct(.keep_all = T)
df.gtop <- dat.gtop %>% select(locus_id,cs,Tissue) %>% distinct(.keep_all = T)
df.cran <- dat.cran %>% select(Gene,CS_ID,Tissue) %>% distinct(.keep_all = T)


df_count.cran <- df.cran %>% group_by(Gene,Tissue) %>% summarise(cs_count=n())
df_count.gtex <- df.gtex %>% group_by(locus_id,Tissue) %>% summarise(cs_count=n())
df_count.gtop <- df.gtop %>% group_by(locus_id,Tissue) %>% summarise(cs_count=n())

df_count.cran$group <- "GTEx+GTOP"
df_count.gtex$group <- "GTEx"
df_count.gtop$group <- "GTOP"

names(df_count.gtex) <- c("Gene","Tissue","cs_count","group")
names(df_count.gtop) <- c("Gene","Tissue","cs_count","group")
df_count <- rbind(df_count.cran,df_count.gtex,df_count.gtop)

# extract genes with single cs
df_count.cran$tissue_gene <- paste(df_count.cran$Gene,df_count.cran$Tissue,sep=":")
df_count.gtex$tissue_gene <- paste(df_count.gtex$Gene,df_count.gtex$Tissue,sep = ":")
df_count.gtop$tissue_gene <- paste(df_count.gtop$Gene,df_count.gtop$Tissue,sep = ":")
x.cran <- as.character(df_count.cran$tissue_gene[df_count.cran$cs_count==1])
x.gtex <- as.character(df_count.gtex$tissue_gene[df_count.gtex$cs_count==1])
x.gtop <- as.character(df_count.gtop$tissue_gene[df_count.gtop$cs_count==1])

overlap_tissuegenes <- intersect(intersect(x.gtex,x.gtop),x.cran)

df.cran <- dat.cran %>% filter(tissue_gene %in% overlap_tissuegenes) %>% select(Gene,CS_ID,MAX_PIP,CS_LENGTH,Tissue,tissue_gene)
df.gtex <- dat.gtex %>% filter(tissue_gene %in% overlap_tissuegenes) %>% select(locus_id,cs,pip,cs_size,Tissue,tissue_gene) %>% 
  group_by(locus_id,cs,Tissue) %>% mutate(max_pip = max(pip)) %>% ungroup() %>% filter(pip==max_pip) %>% 
  select(locus_id,cs,max_pip,cs_size,Tissue,tissue_gene) %>% distinct(.keep_all = T)
df.gtop <- dat.gtop %>% filter(tissue_gene %in% overlap_tissuegenes) %>% select(locus_id,cs,pip,cs_size,Tissue,tissue_gene) %>% 
  group_by(locus_id,cs,Tissue) %>% mutate(max_pip = max(pip)) %>% ungroup() %>% filter(pip==max_pip) %>% 
  select(locus_id,cs,max_pip,cs_size,Tissue,tissue_gene) %>% distinct(.keep_all = T)

df.cran$group <- "GTEx+GTOP"
df.gtex$group <- "GTEx"
df.gtop$group <- "GTOP"
names(df.cran) <- c("locus_id","cs","max_pip","cs_size","Tissue","tissue_gene","group")
df_plot1 <- rbind(df.cran,df.gtex,df.gtop)
df_plot1$group <- factor(df_plot1$group,levels = c("GTOP","GTEx","GTEx+GTOP"))

p2 <- ggplot(df_plot1,aes(x=Tissue,y=log2(cs_size),fill=group)) + 
  geom_boxplot(width=.5) + theme_pubr() + 
  scale_fill_manual(breaks =c("GTOP","GTEx","GTEx+GTOP"), values = c("#B65844","#E2C396","#7784A3") );p2

# Extended Data Fig.7d correlation of PIP between single population fine-mapping and cross-ancestry fine-mapping
x <- df.cran %>% mutate(maxPIP_cran=max_pip,CS_size_cran=cs_size) %>% select(tissue_gene,maxPIP_cran,CS_size_cran)
y <- df.gtex %>% mutate(maxPIP_gtex=max_pip,CS_size_gtex=cs_size) %>% select(tissue_gene,maxPIP_gtex,CS_size_gtex) 
z <- df.gtop %>% mutate(maxPIP_gtop=max_pip,CS_size_gtop=cs_size) %>% select(tissue_gene,maxPIP_gtop,CS_size_gtop) 
xy <- merge(x,y,by="tissue_gene")
xyz <- merge(xy,z,by="tissue_gene")

xyz$maxPIP_gtex_gtop <- apply(xyz[,c(4,6)],1,max)
xyz$maxPIP_all <- apply(xyz[,c(2,4,6)],1,max)

xyz.f <- xyz %>% filter(CS_size_cran>1)
p3.1 <- ggplot(xyz,aes(x=maxPIP_cran,y=maxPIP_gtex)) + geom_point(aes(color=maxPIP_all)) + theme_pubr() + stat_density2d(color="grey")+scale_color_distiller(palette = "OrRd");p3.1
p3.2 <- ggplot(xyz,aes(x=maxPIP_cran,y=maxPIP_gtop)) + geom_point(aes(color=maxPIP_all)) + theme_pubr()+ stat_density2d(color="grey")+scale_color_distiller(palette = "OrRd");p3.2
p3.3 <- ggplot(xyz,aes(x=maxPIP_cran,y=maxPIP_gtex_gtop)) + geom_point(aes(color=maxPIP_all)) + theme_pubr()+ stat_density2d(color="grey") + scale_color_distiller(palette = "OrRd");p3.3
#cowplot::plot_grid(p3.1,p3.2,p3.3,ncol=3,align = "h")

#pdf(file="ext7d.pdf",width = 4.5,height = 6.0)
print(p3.3)
#dev.off()
