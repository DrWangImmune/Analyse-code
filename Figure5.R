#Figure 5
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
#Figure 5A
scobj <- readRDS(file = 'output/GSE157783_micro_annotated.rds')
Idents(scobj) <- 'celltype'
metadata <- scobj@meta.data
dd2 <- data.frame(table(metadata$celltype))
colnames(dd2) <- c("celltype","num")
dd2$Celltype <- paste0(dd2$celltype,paste0(" (",dd2$num,")"))
labels <- dd2$Celltype
color <- brewer.pal(12,'Set1')
DimPlot(scobj, reduction = 'tsne', label = F, label.size = 4.5) +  
  scale_color_manual(values = color,name = 'Celltype',labels = labels) +  
  theme(panel.background = element_blank(),  
        legend.title = element_text(size = 14,color = 'black'),
        legend.text = element_text(size = 12,color = 'black'),
        plot.title = element_text(size = 16,colour = 'black',hjust = 0.5),
        axis.line = element_blank(),  
        axis.title = element_blank(),  
        axis.ticks = element_blank(),  
        axis.text = element_blank()) +  
  labs(title = 'Microglia subtype')+
  annotate("segment", 
           x = -25, xend = -16, 
           y = -35, yend = -35, 
           arrow = arrow(type = "closed", length = unit(0.1, "inches")), 
           color = "black", size = 0.75) +
  annotate("segment",  
           x = -25, xend = -25, 
           y = -35, yend = -25, 
           arrow = arrow(type = "closed", length = unit(0.1, "inches")), 
           color = "black", size = 0.75) + 
  annotate("text",  
           x = -21, y = -37, 
           label = "tSNE1", size = 4, 
           color = "black") +
  annotate("text", 
           x = -27, y = -30, 
           label = "tSNE2", size = 4, 
           color = "black", angle = 90) 
#Figure 5E
markers <- FindAllMarkers(scobj, only.pos = TRUE)
top <- markers |>
  dplyr::group_by(cluster) |>
  dplyr::filter(avg_log2FC > 1) |>
  dplyr::filter(p_val<0.05) |>
  dplyr::slice_head(n = 50) |>
  dplyr::ungroup()
gg <- bitr(top$gene, 'SYMBOL', 'ENTREZID', 'org.Hs.eg.db') 
top <- merge(top, gg, by.x='gene', by.y = 'SYMBOL')
kk <- compareCluster(ENTREZID~cluster, 
                     data = top, 
                     fun="enrichKEGG",
                     organism = 'hsa',
                     pvalueCutoff=0.05)
dotplot(kk, label_format=100,showCategory = 6) + 
  aes(x=sub("\n.*", "", Cluster)) + 
  scale_fill_distiller(palette = "Spectral", direction = 1) +
  labs(x = NULL,title = 'Microglia subtype')+
  theme(axis.text.x = element_text(angle=45, hjust=1),
        axis.text = element_text(size = 12,color = 'black'),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size = 0.8,linetype="solid"),
        plot.title = element_text(size = 16,hjust = 0.5,colour = 'black'),
        legend.frame = element_rect(colour = "black"),          
        legend.ticks = element_line(colour = "black", linewidth  = 0),          
        legend.key.width = unit(0.3, "cm"),          
        legend.key.height = unit(0.6, "cm"),          
        axis.line = element_blank())+
  guides(fill = guide_colorbar(order = 1),shape = guide_legend(order = 2))+
  scale_y_discrete(labels=function(y) str_wrap(y, width=50))
#Figure 5F
#DAM/HM
VlnPlot(scobj, features = c('TMEM119',"P2RY12",'CX3CR1','CSF1R'),cols=brewer.pal(4,'Set3'),ncol = 2)
VlnPlot(scobj, features = c('APOE',"SPP1",'CSTB','LGALS3','TREM2',"CLEC7A"),cols=brewer.pal(4,'Set3'),ncol = 3)
plot_density(scobj,features = c('TMEM119',"P2RY12",'CX3CR1','CSF1R'),reduction = 'tsne')+plot_layout(ncol = 2)
plot_density(scobj,features = c('APOE',"SPP1",'CSTB','LGALS3','TREM2',"CLEC7A"),reduction = 'tsne')+plot_layout(ncol = 3)
#M1/M2
VlnPlot(scobj, features = c('IL10','MRC1',"ARG1","TGFB1"),cols=brewer.pal(6,'Set3'),ncol = 2)
VlnPlot(scobj, features = c('CD68','IL1B','IL6',"CD86"),cols=brewer.pal(6,'Set3'),ncol = 2)
plot_density(scobj,features = c('IL10','MRC1',"ARG1","TGFB1"),reduction = 'tsne')+plot_layout(ncol = 2)
plot_density(scobj,features = c('CD68','IL1B','IL6',"CD86"),reduction = 'tsne')+plot_layout(ncol = 2)