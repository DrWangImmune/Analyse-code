#Figure 6
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
#Figure 6A-D
rm(list = ls())
scobj <- readRDS(file = 'output/micro_annotated.rds')
expr_matrix = GetAssayData(scobj, slot = "counts")
p_data <- scobj@meta.data
p_data$celltype <- scobj@active.ident
f_data <- data.frame(gene_short_name = row.names(scobj),row.names = row.names(scobj))
pd <- new("AnnotatedDataFrame",data = p_data)
fd <- new("AnnotatedDataFrame",data = f_data)
cds <- newCellDataSet(expr_matrix,phenoData = pd,featureData = fd,lowerDetectionLimit = 0.5,expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds,min_expr = 0.1)
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table,mean_expression >= 0.1 & dispersion_empirical >=1*dispersion_fit)$gene_id
cds <- setOrderingFilter(cds,disp.genes)
plot_ordering_genes(cds)
diff <- differentialGeneTest(cds[disp.genes,],fullModelFormulaStr = "~subcelltype",cores = 10)
deg <- subset(diff,qval < 0.01)
deg <- deg[order(deg$qval,decreasing = F),]
ordergene <- rownames(deg)
cds <- setOrderingFilter(cds,ordergene)
ordergene <- row.names(deg)[order(deg$qval)]
cds <- reduceDimension(cds,max_components = 2,method = "DDRTree")
cds <- orderCells(cds)
cds <- orderCells(cds,root_state = GM_state(cds))
plot_cell_trajectory(cds,color_by = "Pseudotime",
                     size = 1,show_branch_points = T,
                     cell_size = 2.5,cell_link_size = 1,
                     show_backbone = T)+ 
  scale_color_distiller(palette = "YlGnBu",direction = 1)+
  labs(x = "Component 1", y = "Component 2", title = 'Pseudotime') + 
  theme(panel.background = element_blank(),  
        panel.border = element_blank(), 
        axis.text = element_text(size = 12, colour = 'black'),  
        axis.title = element_text(size = 14, color = 'black'),  
        plot.title = element_text(size = 16, hjust = 0.5, color = 'black'),  
        legend.text = element_text(size = 12, color = 'black'), 
        legend.title = element_text(size = 12, color = 'black', vjust = 1), 
        legend.frame = element_rect(colour = "black"), 
        legend.ticks = element_line(colour = "black", linewidth = 0), 
        legend.key.width = unit(0.6, "cm"),  
        legend.key.height = unit(0.3, "cm"),  
        legend.position = 'top')  
color <- brewer.pal(8,'Set1')
plot_cell_trajectory(cds,color_by = "State",size = 1,
                     show_branch_points = T,
                     cell_size = 2.5,cell_link_size = 1,
                     show_backbone = T)+
  scale_color_manual(values = color)+
  labs(x = "Component 1", y = "Component 2", title = 'State') + 
  theme(panel.background = element_blank(),  
        panel.border = element_blank(), 
        axis.line = element_line(size = 0.65, color = 'black'), 
        axis.ticks = element_line(size = 0.65, color = 'black'),  
        axis.text = element_text(size = 12, colour = 'black'),  
        axis.title = element_text(size = 14, color = 'black'),  
        plot.title = element_text(size = 16, hjust = 0.5, color = 'black'),  
        legend.text = element_text(size = 12, color = 'black'), 
        legend.title = element_text(size = 12, color = 'black', vjust = 1), 
        legend.frame = element_rect(colour = "black"), 
        legend.ticks = element_line(colour = "black", linewidth = 0), 
        legend.key.width = unit(0.6, "cm"),  
        legend.key.height = unit(0.3, "cm"),  
        legend.position = 'top')  +
  guides(color = guide_legend(nrow = 1,title=NULL)) 
plot_cell_trajectory(cds,color_by = "celltype",size = 1,
                     show_branch_points = T,
                     cell_size = 2.5,cell_link_size = 1,
                     show_backbone = T)+
  scale_color_manual(values = color)+
  labs(x = "Component 1", y = "Component 2", title = "Subcelltype") + 
  theme(panel.background = element_blank(),  
        panel.border = element_blank(), 
        axis.line = element_line(size = 0.65, color = 'black'), 
        axis.ticks = element_line(size = 0.65, color = 'black'),  
        axis.text = element_text(size = 12, colour = 'black'),  
        axis.title = element_text(size = 14, color = 'black'),  
        plot.title = element_text(size = 16, hjust = 0.5, color = 'black'),  
        legend.text = element_text(size = 12, color = 'black'), 
        legend.title = element_text(size = 12, color = 'black', vjust = 1), 
        legend.frame = element_rect(colour = "black"), 
        legend.ticks = element_line(colour = "black", linewidth = 0), 
        legend.key.width = unit(0.6, "cm"),  
        legend.key.height = unit(0.3, "cm"),  
        legend.position = 'top')  +
  guides(color = guide_legend(nrow = 1,title=NULL)) 
color <- brewer.pal(8,'Accent')
plot_cell_trajectory(cds,color_by = "group",size = 1,
                     show_branch_points = T,
                     cell_size = 2.5,cell_link_size = 1,
                     show_backbone = T)+
  scale_color_manual(values = color)+
  labs(x = "Component 1", y = "Component 2", title = 'Group') + 
  theme(panel.background = element_blank(),  
        panel.border = element_blank(), 
        axis.line = element_line(size = 0.65, color = 'black'), 
        axis.ticks = element_line(size = 0.65, color = 'black'),  
        axis.text = element_text(size = 12, colour = 'black'),  
        axis.title = element_text(size = 14, color = 'black'),  
        plot.title = element_text(size = 16, hjust = 0.5, color = 'black'),  
        legend.text = element_text(size = 12, color = 'black'), 
        legend.title = element_text(size = 12, color = 'black', vjust = 1), 
        legend.frame = element_rect(colour = "black"), 
        legend.ticks = element_line(colour = "black", linewidth = 0), 
        legend.key.width = unit(0.6, "cm"),  
        legend.key.height = unit(0.3, "cm"),  
        legend.position = 'top')  +
  guides(color = guide_legend(nrow = 1,title=NULL)) 
#Figure 6E
res <- BEAM(cds[ordergene,],branch_point = 1,cores = 4,progenitor_method = 'duplicate')
res <- res[order(res$qval),]
plot_genes_branched_heatmap(cds[row.names(subset(res,qval < 1e-10)),],
                            branch_point = 1,num_clusters = 3,cores = 4,use_gene_short_name = T,show_rownames = T)
df <- plot_genes_branched_heatmap2(cds[row.names(subset(res,qval < 1e-10)),],
                                   branch_point = 1,num_clusters = 3,cores = 4,
                                   use_gene_short_name = T,show_rownames = T)
enrich <- enrichCluster(object = df,OrgDb = "org.Hs.eg.db",
                        type = "KEGG",organism = "hsa",
                        pvalueCutoff = 0.5,topn = 6,seed = 123)
markGenes = sample(unique(df$wide.res$gene),25,replace = F)
visCluster(object = df,
           plot.type = "both",
           line.col = 'black',
           ht.col.list = list(col_range = seq(-2,2,length.out = 100),
                              col_color = colorRampPalette(c('skyblue',"white",'pink'))(100)),
           pseudotime_col = brewer.pal(3,'Paired'),
           ctAnno.col =brewer.pal(3,'Set2'),
           column_names_rot = 45,show_row_dend = F,
           markGenes = markGenes,markGenes.side = "left",
           annoTerm.data = enrich,
           go.col = rep(brewer.pal(3,'Set2'),each = 6),
           line.side = "left")