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

#用 ggplot 绘图
ggplot(umap1, aes(x = UMAP_1, y = UMAP_2, color = tx), pt.size=0.5) + 
  geom_point(size = 0.2, alpha = 0.6) + 
  scale_color_manual(values = getPalette(colourCount))+
  theme_void()
ggsave(file = "umap.by_cluster.pdf", width = 20, height = 20)
ggsave(file = "umap.by_cluster.png", width = 15, height = 15, units = "in", dpi = 300)

#按celltype绘图并着色
celltype.anno <- read.delim("CellType.txt")
colnames(celltype.anno)[1] = "seurat_clusters"
celltype.anno$seurat_clusters <- as.factor(celltype.anno$seurat_clusters)
c <- data@meta.data %>%
  as.data.frame()%>%
  left_join(celltype.anno)
head(c,3)

umap2 = data@reductions$umap@cell.embeddings %>%
      as.data.frame() %>% cbind(tx = c$`Celltype`)

#颜色梯度
colourCount = length(unique(umap2$tx))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

#用 ggplot 绘图
ggplot(umap2, aes(x = UMAP_1, y = UMAP_2, color = tx)) + 
  geom_point(size = 0.2, alpha = 0.6) + 
  scale_color_manual(values = getPalette(colourCount))+
  theme_void()
ggsave(file = "umap.by_celltype.pdf", width = 20, height = 20)
ggsave(file = "umap.by_celltype.png", width = 15, height = 15, units = "in", dpi = 300)
```
#按组染色

#substring-截取字符
umap3 = data@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% cbind(tx = data@meta.data$orig.ident) %>% 
  mutate(group = substring(tx,1,2))

#用 ggplot 绘图
p <- ggplot(umap3, aes(x = UMAP_1, y = UMAP_2, color = group), pt.size=0.5) + 
  geom_point(size = 0.2, alpha = 0.6) + 
  scale_color_npg()+
  theme_void()
ggsave(file = "umap.by_group.pdf", width = 20, height = 20)
ggsave(file = "umap.by_group.png", width = 20, height = 20, units = "in", dpi = 300)

p4 <- p + facet_wrap( ~ group)
ggsave(file = "umap.by_each_group.pdf", p4, width = 20, height = 6.5)
ggsave(file = "umap.by_each_group.png", p4, width = 30, height = 10, units = "in", dpi = 300)