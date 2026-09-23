#==============================================#
# SMR #
# Extended-Figure-9#
#==============================================#

setwd("/media/london_A/mengxin/GTOP_code/extend/extend_9/")

library(tidygraph)
library(ggraph)

# Extended.Fig.9 SMR and colocalization ------------------------------------

rename<-fread("input/rename_traits_v2.txt")


coloc_res<-fread("input/snv_eqtl_coloc_clean.txt")
coloc_res[which(coloc_res$gwas_name=="593"),]$gwas_name<-"Hematuria"
gsub("@","-",coloc_res$gwas_name)->coloc_res$gwas_name
gsub("-Sakaue-2021","",coloc_res$gwas_name)->coloc_res$gwas_name
gsub("-Shirai-2021","",coloc_res$gwas_name)->coloc_res$gwas_name
coloc_res_e<-coloc_res
coloc_res_e%<>%left_join(rename%>%select(gwas_id,combined),by=c("gwas_name"="gwas_id") )
coloc_res_e%<>%mutate(IID=paste0(combined,"_",xqtl_type,"_",symbol))


SMR_res<-fread("input/snv_eqtl_smr_clean.txt")
gsub("@","-",SMR_res$gwas_name)->SMR_res$gwas_name
gsub("-Sakaue-2021","",SMR_res$gwas_name)->SMR_res$gwas_name
gsub("-Shirai-2021","",SMR_res$gwas_name)->SMR_res$gwas_name
SMR_res_e<-SMR_res%>%filter(xqtl_type=="snv_eqtl")
SMR_res_e%<>%left_join(rename%>%select(gwas_id,combined),by=c("gwas_name"="gwas_id"))
SMR_res_e%<>%mutate(IID=paste0(combined,"_",xqtl_type,"_",symbol))

intersect(unique(coloc_res_e$gwas_name),unique(SMR_res_e$gwas_name)) #119 traits
both_coloc_SMR<-intersect(coloc_res_e$IID,SMR_res_e$IID)

coloc_res_e%>%filter(IID%in%both_coloc_SMR)%>%dplyr::select(gwas_name,disease_trait,trait_class,symbol,combined)%>%unique()%>%mutate(color="#f27830")->coloc_SMR_sim
coloc_res_e%>%filter(!IID%in%both_coloc_SMR)%>%dplyr::select(gwas_name,disease_trait,trait_class,symbol,combined)%>%unique()%>%mutate(color="#1d73b6")->coloc_only_sim
SMR_res_e%>%filter(!IID%in%both_coloc_SMR)%>%dplyr::select(gwas_name,disease_trait,trait_class,symbol,combined)%>%unique()%>%mutate(color="#6dc2d2")->SMR_only_sim


rbind(coloc_SMR_sim,coloc_only_sim,SMR_only_sim)->overall_sim
overall_sim$combined<-gsub("_"," ",overall_sim$combined)
overall_sim%<>%mutate(node=paste0(trait_class,"/",combined,"/",symbol))

overall_sim$combined<-gsub("_"," ",overall_sim$combined)
data_input<-overall_sim%>%dplyr::select(Category=trait_class,trait=combined,Gene=symbol)%>%unique()%>%filter(Category!="Quantitative_trait")

source("multi-circle_function.R")
data_input$ID <- 1

index_level <- c("Category","trait","Gene")
nodes_data_input <- gather_graph_node(data_input,index=index_level,root="eQTL")
nodes_data_input$trait <- nodes_data_input$node.branch
head(nodes_data_input)

edges_data_input <- gather_graph_edge(data_input,index=index_level, root="eQTL")
edges_data_input$trait <- stringr::str_split_i(edges_data_input$from, "/",1)
head(edges_data_input)

graph_data_input <- tbl_graph(nodes_data_input,edges_data_input)

# my_col<-c(
#   "quantitative_trait" = "#2F4F4F",
#   "endocrine_metabolic" = "#4682B4",
#   "immune" = "#9e4832",
#   "circulatory_system" = "#487c51",
#   "neoplasms" = "#7371a3",
#   "dermatologic" = "#7690a4",
#   "genitourinary" = "#B0C4DE",
#   "digestive" = "#f0d795",
#   "hematopoietic" = "#00008B",
#   "sensory" = "#c4598a",
#   "nervous_system" = "#7c776b",
#   "respiratory" = "#5192c1",
#   "musculoskeletal" = "#DB7093",
#   "symptoms" = "#2F4F4F"
# )

my_col<-c(c(
  "Quantitative_trait" = "#2F4F4F",
  "Circulatory_system" = "#487c51",
  "Dermatologic" = "#7690a4",
  "Digestive" = "#f0d795",
  "Endocrine_metabolic" = "#4682B4",
  "Genitourinary" = "#B0C4DE",
  "Infectious_diseases" = "#9e4832",
  "Mental_disorders" = "#7c776b",
  "Musculoskeletal" = "#DB7093",
  "Neoplasms" = "#7371a3",
  "Neurological" = "#00008B",
  "Respiratory" = "#5192c1",
  "Sense_organs" = "#c4598a",
  "Symptoms" = "#2F4F4F",
  "eQTL" = "#9ac294",
  "sQTL" = "#7b86a7",
  "eQTL,sQTL" = "#bc9cc1", Hematopoietic ="#00bfc4"
))

gm <- ggraph(graph_data_input,layout = 'dendrogram', circular = TRUE)+
  geom_edge_diagonal(aes(color=trait),
                     alpha=1/4,
                     show.legend = F) +
  geom_node_point(aes(size=node.size,
                      color=trait),
                  alpha=1/3,
                  show.legend = T,
  ) +
  guides(size ="none") + 
  coord_fixed() + 
  theme(legend.position = "right",
        panel.background = element_rect(fill = NA, color = NA),
        plot.background = element_rect(fill = NA, color = NA))+
  scale_size(range = c(1,10))+
  scale_color_manual(values = my_col)+        
  ggraph::scale_edge_color_manual(values = my_col)

gm


gene_dat <- filter(gm$data, node.level=="Gene")
gene_dat$node <- gene_dat$node.name
gene_dat%<>%left_join(overall_sim%>%dplyr::select(node,color)%>%unique(),by="node")

#gene_dat$color<-"black"
gm1<-gm +
  geom_node_text(
    aes(
      x = 1.0175 * x, 
      y = 1.0175 * y,
      label = node.short_name,
      angle = -((-node_angle(x, y) + 90) %% 180) + 90),
    color = gene_dat$color,
    fontface="bold",
    size = 1.55,
    hjust = 'outward',
    data =gene_dat)
gm1

trait_dat <- filter(gm$data, node.level=="trait")
trait_dat$node.short_name<-gsub("_"," ",trait_dat$node.short_name)

gm2 <- gm1+
  geom_node_text(
    aes(x = 1.0175 * x,
        y = 1.0175 * y,
        label=node.short_name,
        angle = -((-node_angle(x, y) + 90) %% 180) + 90),
    color="black",
    #fontface="bold",
    size=1.5,
    hjust = 'outward',
    vjust=-0.8,
    family="sans", # 等价于Arial
    data = trait_dat
  )

gm2

gm3 <- gm2 +
  geom_node_text(
    aes(label=gsub("_"," ",node.short_name),
        angle = -((-node_angle(x, y) + 90) %% 180) + 90,
        color=trait),
    fontface="bold",
    show.legend = F,
    size=2.5,
    family="sans",
    data = filter(gm$data, node.level=="Category")
  )

gm3



