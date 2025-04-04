#Figure 1
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
#Figure 1A
scobj <- readRDS(file = 'output/GSE157783_annotated.rds')
Idents(scobj) <- 'celltype'
metadata <- scobj@meta.data
dd2 <- data.frame(table(metadata$celltype))
colnames(dd2) <- c("celltype","num")
dd2$Celltype <- paste0(dd2$celltype,paste0(" (",dd2$num,")"))
labels <- dd2$Celltype
color <- brewer.pal(12,'Paired')
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
  labs(title = 'GSE157783 Midbrain')+
  annotate("segment", 
           x = -40, xend = -30, 
           y = -40, yend = -40, 
           arrow = arrow(type = "closed", length = unit(0.1, "inches")), 
           color = "black", size = 0.75) +
  annotate("segment",  
           x = -40, xend = -40, 
           y = -40, yend = -30, 
           arrow = arrow(type = "closed", length = unit(0.1, "inches")), 
           color = "black", size = 0.75) + 
  annotate("text",  
           x = -35, y = -42.5, 
           label = "tSNE1", size = 4, 
           color = "black") +
  annotate("text", 
           x = -42.5, y = -36, 
           label = "tSNE2", size = 4, 
           color = "black", angle = 90) 
#Figure 1B
data <- as.data.frame(table(scobj$group,scobj$celltype))
colnames(data) <- c("group","CellType","Freq")
data <- data %>% 
  group_by(group) %>% 
  mutate(Total = sum(Freq)) %>% 
  ungroup() %>% 
  mutate(Percent = Freq/Total) %>% 
  as.data.frame()
data$group <- factor(data$group,levels =  c("Control","PD"))
ggplot(data, aes(x = group, y = Percent, fill = CellType)) +
  scale_fill_manual(values = color)+
  geom_bar(position = "fill", stat="identity", 
           color = 'white', alpha = 0.9, width = 0.8) +
  scale_y_continuous(labels = scales::percent_format(),expand = c(0, 0))+
  labs(title = 'Midbrain cellratio',x = NULL,y = NULL)+
  theme(panel.background = element_rect(fill = 'transparent'), 
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 16,hjust = 0,face = 'bold',vjust = 2,colour = 'black'),
        legend.text = element_text(size = 14,color = "black"),
        axis.text.y = element_text(size = 14,color = "black"),
        axis.text.x = element_text(size = 14,color = "black",angle = 45,hjust = 1)) + 
  guides(fill=guide_legend(title=NULL))
#Figure 1C
markers.all<-FindAllMarkers(scobj,only.pos=TRUE)
top <- markers.all %>% group_by(cluster) %>% top_n(2, avg_log2FC)
genes<- c('PLP1',"MOBP",'SLC17A6','FSTL4','FCGR1A','CD74',
          'PDGFRA',"VCAN","AQP4","GFAP","CLDN5","PECAM1",
          "COL1A2","PDGFRB","GATA3","GAD1","FOXJ1","ADGB",
          "FGF5","FAT2","GRIK1","DMBX1","TH","SLC6A3")
aver <- AverageExpression(scobj,features= genes,group.by = 'celltype',slot= 'data')
aver <- as.data.frame(aver$RNA)
marker<- data.frame(marker = top$cluster,row.names = genes)
celltype<- data.frame(celltype = colnames(aver), row.names = colnames(aver))
names(color) <- celltype$celltype
anno_col<- list(celltype = color,marker= color)
pheatmap(as.matrix(aver),
         scale= 'row',
         show_colnames = F,
         cluster_rows= FALSE,
         cluster_cols= FALSE,
         annotation_names_row = FALSE,
         annotation_names_col = FALSE,
         annotation_col= celltype,
         annotation_row= marker,
         annotation_colors= anno_col, 
         color= colorRampPalette(c("#1F78B4","white", "#E31A1C"))(100),
         border_color= 'white',fontsize = 11) 
#Figure 1D
rm(list = ls())
scobj <- readRDS(file = 'output/GSE243639_annotated.rds')
Idents(scobj) <- 'celltype'
color <- brewer.pal(12,'Paired')
metadata <- scobj@meta.data
dd2 <- data.frame(table(metadata$celltype))
colnames(dd2) <- c("celltype","num")
dd2$Celltype <- paste0(dd2$celltype,paste0(" (",dd2$num,")"))
labels <- dd2$Celltype
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
  labs(title = 'GSE243639 Substantia nigra compacta')+
  annotate("segment", 
           x = -40, xend = -30, 
           y = -40, yend = -40, 
           arrow = arrow(type = "closed", length = unit(0.1, "inches")), 
           color = "black", size = 0.75) +
  annotate("segment",  
           x = -40, xend = -40, 
           y = -40, yend = -30, 
           arrow = arrow(type = "closed", length = unit(0.1, "inches")), 
           color = "black", size = 0.75) + 
  annotate("text",  
           x = -35, y = -42.5, 
           label = "tSNE1", size = 4, 
           color = "black") +
  annotate("text", 
           x = -42.5, y = -36, 
           label = "tSNE2", size = 4, 
           color = "black", angle = 90) 
#Figure 1E
data <- as.data.frame(table(scobj$group,scobj$celltype))
colnames(data) <- c("group","CellType","Freq")
data <- data %>% 
  group_by(group) %>% 
  mutate(Total = sum(Freq)) %>% 
  ungroup() %>% 
  mutate(Percent = Freq/Total) %>% 
  as.data.frame()
data$group <- factor(data$group,levels =  c("Control","PD"))
ggplot(data, aes(x = group, y = Percent, fill = CellType)) +
  scale_fill_manual(values = color)+
  geom_bar(position = "fill", stat="identity", color = 'white', alpha = 0.9, width = 0.8) +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = c(0, 0))+
  labs(title = 'Substantia nigra compacta cellratio',x = NULL,y = NULL)+
  theme(panel.background = element_rect(fill = 'transparent'), 
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 14,hjust = 0,face = 'bold',
                                  vjust = 2,colour = 'black'),
        legend.text = element_text(size = 14,color = "black"),
        axis.text.y = element_text(size = 14,color = "black"),
        axis.text.x = element_text(size = 14,color = "black",angle = 45,hjust = 1)) + 
  guides(fill=guide_legend(title=NULL))
#Figure 1F
markers.all<-FindAllMarkers(scobj,only.pos=TRUE)
top <- markers.all %>% group_by(cluster) %>% top_n(3, avg_log2FC)
color <- brewer.pal(7,'Paired')
genes<- c("AQP4","GFAP","SLC1A3",'PLP1',"MOBP","ST18",
          'FCGR1A','CD74',"CD163","SNAP25","SYT1","NRG1",
          "DCN","FLT1","CLDN5","TRAC","CD3E","IL7R",
          'PDGFRA',"VCAN","GPR17")
aver <- AverageExpression(scobj,features= genes, group.by = 'celltype',slot= 'data')
aver <- as.data.frame(aver$RNA)
marker<- data.frame(marker = top$cluster,row.names = genes)
celltype<- data.frame(celltype = colnames(aver),row.names = colnames(aver))
names(color) <- celltype$celltype
anno_col<- list(celltype = color,marker= color)
pheatmap(as.matrix(aver),
         scale= 'row',
         show_colnames = F,
         cluster_rows= FALSE,
         cluster_cols= FALSE,
         annotation_names_row = FALSE,
         annotation_names_col = FALSE,
         annotation_col= celltype,
         annotation_row= marker,
         annotation_colors= anno_col, 
         color= colorRampPalette(c("#1F78B4","white", "#E31A1C"))(100),
         border_color= 'white',fontsize = 11) 