#Figure 4
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
#Figure 4C AddModuleScore
rm(list = ls())
scobj <- readRDS(file = 'output/GSE157783_annotated.rds')
genesets <- read.gmt("resource/h.all.v2023.2.Hs.symbols.gmt") 
inflam <- subset(geneset, term=="INFLAMMATORY RESPONSE")
inflam <- list(inflam$gene)
names(inflam)[1] <- 'Inflammatory Response'
DefaultAssay(scobj) <- 'RNA'
scobj <- AddModuleScore(scobj,features=inflam,name = "score",assay = "RNA")
metadata <- scobj@meta.data
ggplot(scobj@meta.data, aes(x=fct_reorder(celltype,score1,.desc = T), y=score1, fill=celltype)) + 
  scale_fill_manual(values = brewer.pal(12,'Paired'))+
  geom_boxplot(linewidth = 0.8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        plot.title = element_text(size = 16,hjust = 0.5,colour = 'black'),
        panel.background = element_blank(),
        axis.line = element_line(color = 'black',size = 0.6),
        axis.ticks = element_line(color = 'black',size = 0.6),
        axis.text = element_text(color = 'black',size = 12),
        axis.title.y = element_text(color = 'black',size = 16)) + 
  xlab(NULL) + 
  labs(title = 'Midbrain AddModuleScore')+
  ylab("Inflammatory Response Score")
#Figure 4E UCell
rm(list = ls())
scobj <- readRDS(file = 'output/GSE157783_annotated.rds')
geneset <- read.gmt("resource/h.all.v2023.2.Hs.symbols.gmt") 
unique(geneset$term)
inflam <- subset(geneset, term=="INFLAMMATORY RESPONSE")
inflam <- list(inflam$gene)
names(inflam)[1] <- 'Inflammatory Response'
DefaultAssay(scobj) <- 'RNA'
scobj <- AddModuleScore_UCell(scobj,features=inflam,name="_inflam_score")
original_name = "Inflammatory Response_inflam_score"
expected_name = "score1"
location = grep(original_name,colnames(scobj@meta.data))
colnames(scobj@meta.data)[location]= expected_name
ggplot(scobj@meta.data, aes(x=fct_reorder(celltype,score1,.desc = T), y=score1, fill=celltype)) + 
  scale_fill_manual(values = brewer.pal(12,'Paired'))+
  geom_boxplot(linewidth = 0.8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        plot.title = element_text(size = 16,hjust = 0.5,colour = 'black'),
        panel.background = element_blank(),
        axis.line = element_line(color = 'black',size = 0.6),
        axis.ticks = element_line(color = 'black',size = 0.6),
        axis.text = element_text(color = 'black',size = 12),
        axis.title.y = element_text(color = 'black',size = 16)) + 
  xlab(NULL) + 
  labs(title = 'Midbrain UCell Score')+
  ylab("Inflammatory Response Score")
#Figure 4G AUCell
rm(list = ls())
scobj <- readRDS(file = 'output/GSE157783_annotated.rds')
geneset <- read.gmt('resource/h.all.v2023.2.Hs.symbols_modify.gmt')
geneset <- split(geneset$gene,geneset$term)
geneset <- geneset[c("INFLAMMATORY RESPONSE")]
exprMatrix <- GetAssayData(scobj,assay = 'RNA',layer = 'data')
plotGeneCount(exprMatrix)
cell_rankings <- AUCell_buildRankings(exprMatrix,splitByBlocks = T)
cell_AUC <- AUCell_calcAUC(geneSets = geneset,cell_rankings,aucMaxRank = nrow(cell_rankings)*0.1)
cell_assigment <- AUCell_exploreThresholds(cell_AUC,plotHist = T,assignCells = T)
cell_assigment$`INFLAMMATORY RESPONSE`$aucThr
warningMag <- sapply(cell_assigment,function(x) x$aucThr$comment)
warningMag[which(warningMag!='')]
sapply(cell_assigment,function(x) x$aucThr$selected)
geneSetName <- rownames(cell_AUC)[grep("INFLAMMATORY RESPONSE", rownames(cell_AUC))]
AUCell_plotHist(cell_AUC[geneSetName,], aucThr=0.09377553)
abline(v=0.09377553)
cellsTsne<-Embeddings(object=scobj,reduction="tsne")
cellsTsne<-cellsTsne[,c(1,2)]
colnames(cellsTsne)<-c("tSNE1","tSNE2")
selectedThresholds <- getThresholdSelected(cell_assigment)
selectedThresholds[1] <-  0.09377553
AUCell_plotTSNE(tSNE=cellsTsne, exprMat=exprMatrix, cellsAUC=cell_AUC[1,], thresholds=selectedThresholds)
geneSet <- "INFLAMMATORY RESPONSE"
aucs <- as.numeric(getAUC(cell_AUC)[geneSet,])
scobj$"INFLAMMATORY RESPONSE AUCell" <- aucs
ggplot(scobj@meta.data, aes(x=fct_reorder(celltype,score2,.desc = T), y=score2, fill=celltype)) + 
  scale_fill_manual(values = brewer.pal(12,'Paired'))+
  geom_boxplot(linewidth = 0.8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        plot.title = element_text(size = 16,hjust = 0.5,colour = 'black'),
        panel.background = element_blank(),
        axis.line = element_line(color = 'black',size = 0.6),
        axis.ticks = element_line(color = 'black',size = 0.6),
        axis.text = element_text(color = 'black',size = 12),
        axis.title.y = element_text(color = 'black',size = 16)) + 
  xlab(NULL) + 
  labs(title = 'Midbrain AUCell Score')+
  ylab("Inflammatory Response Score")
