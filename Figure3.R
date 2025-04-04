#Figure 3
rm(list = ls())
#Load Packages
library(Seurat)
library(harmony)
library(clustree)
library(UCell)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggridges)
library(ggrepel)
library(gghalves)
library(gtools)
library(ggalluvial)
library(pheatmap)
library(ggnetwork)
library(monocle)
library(CellChat)
library(GseaVis)
library(ClusterGVis)
library(Nebulosa)
library(RColorBrewer)
library(tidydr)
library(ggforce)
library(ggrastr)
library(viridis)
library(colorspace)
library(clusterProfiler)
library(cowplot)
library(patchwork)
library(reshape2)
library(GSVA)
library(GSEABase)
#Figure 3A
rm(list = ls())
scobj <- readRDS(file = 'output/GSE157783_annotated.rds')
scdata <- GetAssayData(scobj, slot = "counts")
scdata <- apply(scdata,2,function(xx)xx/sum(xx)*1e4)
scdata <- log2(scdata + 1)
scobj[['RNA']] <- CreateAssayObject(counts = scdata)
rm(list = ls(pattern = 'scdata'))
expr <- AverageExpression(scobj,group.by = 'group',assay = "RNA",slot = 'data')
expr <- expr[[1]]
expr <- expr[rowSums(expr)>0,]
expr <- as.matrix(expr)
cg = names(tail(sort(apply(expr,1,sd)),1000))
pheatmap::pheatmap(cor(expr[cg,]))
geneset <- read.gmt('resource/h.all.v2023.2.Hs.symbols.gmt')
geneset <- split(geneset$gene,geneset$term)
geneset <- lapply(geneset,unique)
geneset <- GeneSetCollection(mapply(function(geneIds,keggId){
  GeneSet(geneIds,geneIdType = EntrezIdentifier(),
          collectionType = KEGGCollection(keggId),
          setName = keggId)
},geneset,names(geneset)))
gsva_result <- gsva(as.matrix(expr),geneset)
gsva <- melt(gsva_result,id.vars = 'Genesets')
colnames(gsva) <- c("hallmark","celltype","Score")
gsva$celltype <- factor(gsva$celltype,levels = c("PD",'Control'))
ggplot(gsva,aes(x = celltype, y = hallmark,color = Score))+
  geom_point(aes(color=Score),size = 5, alpha = 0.8)+
  theme(panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 16,hjust = 0.5,color = 'black',face = 'bold'),
        panel.background = element_rect(color = 'black',size = 1,fill = 'transparent'),
        axis.text.x=element_text(size=9,color="black",angle=45,hjust = 1), 
        axis.text.y=element_text(size=12,color="black",face = 'bold'),
        legend.frame = element_rect(colour = "black"),          
        legend.ticks = element_line(colour = "black", linewidth  = 0),          
        legend.key.width = unit(0.3, "cm"),          
        legend.key.height = unit(0.6, "cm"),          
        legend.title = element_text(color = 'black', face = "bold", size=8))+
  scale_color_gradientn(values = seq(0,1,0.1),colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  labs(x=NULL,y=NULL,title = 'GSE157783 Midbrain')+
  coord_flip()
#Figure 3B
rm(list = ls())
scobj <- readRDS(file = 'output/GSE243639_annotated.rds')
scdata <- GetAssayData(scobj, slot = "counts")
scdata <- apply(scdata,2,function(xx)xx/sum(xx)*1e4)
scdata <- log2(scdata + 1)
scobj[['RNA']] <- CreateAssayObject(counts = scdata)
rm(list = ls(pattern = 'scdata'))
expr <- AverageExpression(scobj,group.by = 'group',assay = "RNA",slot = 'data')
expr <- expr[[1]]
expr <- expr[rowSums(expr)>0,]
expr <- as.matrix(expr)
cg = names(tail(sort(apply(expr,1,sd)),1000))
pheatmap::pheatmap(cor(expr[cg,]))
geneset <- read.gmt('resource/h.all.v2023.2.Hs.symbols.gmt')
geneset <- split(geneset$gene,geneset$term)
geneset <- lapply(geneset,unique)
geneset <- GeneSetCollection(mapply(function(geneIds,keggId){
  GeneSet(geneIds,geneIdType = EntrezIdentifier(),
          collectionType = KEGGCollection(keggId),
          setName = keggId)
},geneset,names(geneset)))
gsva_result <- gsva(as.matrix(expr),geneset)
gsva <- melt(gsva_result,id.vars = 'Genesets')
colnames(gsva) <- c("hallmark","celltype","Score")
gsva$celltype <- factor(gsva$celltype,levels = c("PD",'Control'))
ggplot(gsva,aes(x = celltype, y = hallmark,color = Score))+
  geom_point(aes(color=Score),size = 5, alpha = 0.8)+
  theme(panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 16,hjust = 0.5,color = 'black',face = 'bold'),
        panel.background = element_rect(color = 'black',size = 1,fill = 'transparent'),
        axis.text.x=element_text(size=9,color="black",angle=45,hjust = 1), 
        axis.text.y=element_text(size=12,color="black",face = 'bold'),
        legend.frame = element_rect(colour = "black"),          
        legend.ticks = element_line(colour = "black", linewidth  = 0),          
        legend.key.width = unit(0.3, "cm"),          
        legend.key.height = unit(0.6, "cm"),          
        legend.title = element_text(color = 'black', face = "bold", size=8))+
  scale_color_gradientn(values = seq(0,1,0.1),colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  labs(x=NULL,y=NULL,title = 'GSE243639 Substantia nigra compacta')+
  coord_flip()



























