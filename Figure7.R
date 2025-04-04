#Figure 7
rm(list = ls())
#Load Packages
library(Seurat)
library(harmony)
library(clustree)
library(UCell)
library(AUCell)
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
library(SCopeLoomR)
#Figure 7A
rm(list = ls())
scobj <- readRDS(file = 'output/micro_annotated.rds')
geneset <- read.gmt('output/micro.regulons.txt')
geneset$gene <- as.character(geneset$gene)
for(i in 1:nrow(geneset)){
  print(i)
  term = geneset$gene[i]
  term = unlist(strsplit(term, split=",", fixed=T))
  term = paste(term, collapse=" ")
  geneset$gene[i] = term
}
gmt_data <- geneset
colnames(gmt_data) <- c('term','genes')
save_as_gmt <- function(gmt_data, file_name) {
  con <- file(file_name, "w")
  for (i in 1:nrow(gmt_data)) {
    gene_set_name <- as.character(gmt_data$term[i])
    genes <- as.character(gmt_data$genes[i])
    line <- paste(c(gene_set_name, "Description", genes), collapse = "\t")
    writeLines(line, con)
  }
  close(con)
}
regulons <- clusterProfiler::read.gmt("resource/micro_pyscenic.gmt")
rg.names <- unique(regulons$term)
regulon.list <- lapply(rg.names, function(rg) {
  subset(regulons, term == rg)$gene
})
names(regulon.list) <- sub("[0-9]+g", "\\+", rg.names)
summary(sapply(regulon.list, length))
print(regulon.list[1])
sce_SCENIC<-open_loom("output/SCENIC.loom")
regulons_incidMat<-get_regulons(sce_SCENIC,column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(sce_SCENIC, column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(sce_SCENIC)
scdata <- GetAssayData(scobj, slot = "counts")
scdata <- apply(scdata,2,function(xx)xx/sum(xx)*1e4)
scdata <- log2(scdata + 1)
mc.mat <- scdata
cellinfo <- scobj@meta.data[,c('celltype','group',"nFeature_RNA","nCount_RNA")]
colnames(cellinfo)=c('celltype', 'group','nGene' ,'nUMI')
cellTypes <-  as.data.frame(subset(cellinfo,select = 'celltype'))
selectedResolution <- "celltype"
sub_regulonAUC <- regulonAUC
rss <- calcRSS(AUC=getAUC(sub_regulonAUC),cellAnnotation=cellTypes[colnames(sub_regulonAUC),selectedResolution])
B_rss <- as.data.frame(rss)
celltype <- "Microglia 1"
rssRanklist <- list()
for(i in 1:length(celltype)) {
  data_rank_plot <- cbind(as.data.frame(rownames(B_rss)),
                          as.data.frame(B_rss[,celltype[i]]))
  colnames(data_rank_plot) <- c("TF", "celltype")
  data_rank_plot=na.omit(data_rank_plot)
  data_rank_plot <- data_rank_plot[order(data_rank_plot$celltype,decreasing=T),]
  data_rank_plot$rank <- seq(1, nrow(data_rank_plot))
  
  p <- ggplot(data_rank_plot, aes(x=rank, y=celltype)) + 
    geom_point(size=3, shape=16, color="#1F77B4",alpha =0.4)+
    geom_point(data = data_rank_plot[1:10,],
               size=3, color='#DC050C')+ 
    theme(axis.title = element_text(colour = 'black', size = 14),
          axis.text = element_text(colour = 'black', size = 12),
          axis.line = element_line(colour = 'black',size = 0.6),
          plot.title = element_text(size = 16,hjust = 0.5,color = 'black'),
          panel.background = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    labs(x='Regulons Rank', y='Specificity Score',title =celltype[i])+
    geom_text_repel(data= data_rank_plot[1:10,],
                    aes(label=TF), color="black", size=4, 
                    arrow = arrow(ends="first", length = unit(0.01, "npc")), box.padding = 0.1,
                    segment.size = 0.3, force = 1, max.iter = 3e3)
  rssRanklist[[i]] <- p
}