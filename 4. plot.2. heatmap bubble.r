#加载包
library(Seurat)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
#加载数据
data <- readRDS("macacards/snRNA.2.rds")

#以防出错，先导入，在整合至新表，在将新表cbind过去
celltype.anno <- read.delim("CellType.txt")
colnames(celltype.anno)[1] = "seurat_clusters"
c <- left_join(data@meta.data, celltype.anno)

data@meta.data <- cbind(data@meta.data, c[,7])
colnames(data@meta.data) [7] <- "Cell type"
umap1 = data@reductions$umap@cell.embeddings %>%
      as.data.frame() %>% cbind(tx = data@meta.data$seurat_clusters)

#color
colourCount = length(unique(umap1$tx))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

#gene markers
memory.limit(1000000)
clus.marker <- FindAllMarkers(data)
write.table(clu.marker,"GraphClust.AllMarkerGenes.txt")
clus.marker <- clus.marker %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_logFC)
pdf(file = "heatmap cluster.pdf", width = 9, height = 9)
png(file = "heatmap cluster.png", width = 15, height = 15, units = "in", res = 300)
DoHeatmap(data, features = clus.marker$allmarkers.gene)+ NoLegend()
dev.off()

#cell markers
memory.limit(1000000)
my_levels <- read.delim("culs level.txt", quote= "") $order %>% as.character()
cellmarkers <- read.delim("filtered cell markers.txt", quote= "") $markers %>% as.character()
Idents(data) <- factor(Idents(data), levels= my_levels)
png(file = "heatmap cell1.png", width = 15, height = 8, units = "in", res = 300)
DoHeatmap(data, features = cellmarkers, label = F,slot = "scale.data", draw.lines = F, size = 5.5, hjust = 0) 
ggsave("heatmap cell1.pdf", device = "pdf", width = 23,height = 12, units = "cm" )
dev.off()

#bubble plot，由上往下看，顺序要反过来
my_levels <- read.delim("culs level bubble.txt", quote= "") $order %>% as.character()
Idents(data) <- factor(Idents(data), levels= my_levels)
DotPlot(data, features = cellmarkers, cols = c("#FEE0D2", "#CB181D"), scale.min = 10) + 
  scale_x_discrete("")+
  scale_y_discrete("")+
  RotatedAxis()
ggsave(file = "Bub cluster.png", width = 10, height = 8, limitsize = FALSE)
ggsave(file = "Bub cluster.pdf", width = 10, height = 8, limitsize = FALSE)