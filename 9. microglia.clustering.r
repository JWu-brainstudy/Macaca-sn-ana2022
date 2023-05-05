library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(tidyverse)
library(monocle)
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(openxlsx)
library(ggsci)

data1 <- readRDS("microglia.seuset.rds")


DimHeatmap(data, dims = 1:20, cells = 1000, balanced = TRUE)
JackStrawPlot(data, dims = 1:20)

data1 <- FindNeighbors(data, dims = 1:8)
data1 <- FindClusters(data1, resolution = 0.8)
data1 <- RunUMAP(data1, dims = 1:8)

umap1 = data1@reductions$umap@cell.embeddings %>%
      as.data.frame() %>% cbind(tx = data1@meta.data$seurat_clusters)

#颜色梯度
colourCount = length(unique(umap1$tx))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

ggplot(umap1, aes(x = UMAP_1, y = UMAP_2, color = tx), pt.size=1) + 
  geom_point(size = 1, alpha = 0.6) + 
  scale_color_manual(values = getPalette(colourCount))+
  theme_void()
ggsave(file = "umap.by_cluster.pdf", width = 4, height = 4)
ggsave(file = "umap.by_cluster.png", width = 4, height = 4, units = "in", dpi = 300)

colourCount = length(unique(umap3$tx))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

p <- ggplot(umap3, aes(x = UMAP_1, y = UMAP_2, color = tx), pt.size=1) + 
  geom_point(size = 1, alpha = 0.6) + 
  scale_color_manual(values = getPalette(colourCount))+
  gait.errorcheck()+
  facet_wrap( ~ group)+
  theme_void()
ggsave(file = "umap.by_each_sample.pdf", width = 10, height = 3.5)
ggsave(file = "umap.by_each_sample.png", width = 10, height = 3.5, units = "in", dpi = 300)

#subcluster proportion
cell.info.total <- read.xlsx("subcluster.total.empty.xlsx")
cell.info.total$seurat_clusters<- as.factor(cell.info.total$seurat_clusters)

cell.info <- data1@meta.data %>%
  group_by(orig.ident)%>%
  mutate(all.cells = n()) %>%
  group_by(orig.ident, seurat_clusters, all.cells)%>%
  summarise(cells = n())%>%
  mutate(cell.per = cells/all.cells)%>%
  group_by(orig.ident)%>%
  mutate(group = substring(orig.ident,1,2)) %>% 
  as.data.frame() %>% 
  left_join(cell.info.total,.)

cell.info$cell.per[is.na(cell.info$cell.per)==T] = 0

my.boxpl <- function(df){
  my_comparisons <- list(c("De", "HC"), c("De", "Re"), c("HC", "Re"))
  p = ggboxplot(df, x="group", y="cell.per", color = "group",
               palette = c("#FC4E07","#00AFBB", "#E7B800"), 
               add = "jitter")+
    theme(legend.title=element_blank())+
  stat_compare_means(method = "wilcox.test", comparisons = my_comparisons)+
  stat_compare_means(label.y = 1)
  p
}


p <- ggplot(data = cell.info)


p4 <- p + stat_boxplot(geom = "errorbar", 
                       width=0.4, 
                       aes(x = seurat_clusters, y = cell.per,fill=group), 
                       position = position_dodge(width = 0.6),
                       outlier.shape = NA)

p5 <- p4 + stat_boxplot(geom = "boxplot", 
                       width=0.4, 
                       aes(x = seurat_clusters, y = cell.per,fill=groupz), 
                       position = position_dodge(width = 0.6),
                       outlier.shape = NA)

# p6 <- p5 + geom_point(aes(x = seurat_clusters, y = cell.per,fill=group), width = 0.2)

p1 <- p5 + theme_classic() + theme(legend.position="none")

(p1)                       
                       
```


#5.4 小胶质细胞相关marker
```{r}
# cellmarkers <- read.delim("filtered cell markers.txt", quote= "") $markers %>% as.character()
FeaturePlot(data1, features = "PTPRC", pt.size = 0.5, reduction = "umap")+
  theme_void()
ggsave(file = "markers.ptprc.png", width = 4, height = 4, limitsize = FALSE)
```

#5.5 20230106 小胶质细胞相关marker
```{r}
color <- c("#EAEAEA", "#e41a1c", "#3f007d")
p_iq <- FeaturePlot(data1, features = c("IQGAP2", "FYN", "PDE7A",  "ARHGEF3","BCL11B"), cols = color, pt.size = 0.5, reduction = "umap", ncol = 5)+ theme_void()
ggsave(file = "markers.IQGAP2.png", width = 20, height = 4, limitsize = FALSE)
```



#6. monocle拟时序
#6.1 转换monocle对象
```{r}
x <- as(as.matrix(data1@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = data1@meta.data)
fData <- data.frame(gene_short_name = row.names(x), row.names = row.names(x))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(x,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())


mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)

rm(x,pd,fData,fd)
```

#6.2 seurat marker gene
```{r}
diff.wilcox = FindAllMarkers(data1)
disp.genes <- row.names (subset(diff.wilcox, p_val_adj < 0.01))
mycds <- setOrderingFilter(mycds, disp.genes)
```

#6.3 其他方法计算markergene
```{r}
diff_test_res <- differentialGeneTest(HSMM_myo[expressed_genes,],
              fullModelFormulaStr = "~Media")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)
rm(x)
```


```{r}
#降维
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
#排序
mycds <- orderCells(mycds)
```

按组染色
```{r}
plot4.dep <- plot_cell_trajectory(mycds, color_by = "orig.ident") + 
  scale_color_manual(values = c(rep("#E41A1C",7), rep("NA",7), rep("NA",5)))
plot4.res <- plot_cell_trajectory(mycds, color_by = "orig.ident") + 
  scale_color_manual(values = c(rep("NA",7), rep("NA",7), rep("#4DAF4A",5)))
plot4.HC <- plot_cell_trajectory(mycds, color_by = "orig.ident") + 
  scale_color_manual(values = c(rep("NA",7), rep("#377EB8",7), rep("NA",5)))

plot4 <- plot4.dep|plot4.HC|plot4.res
plot4 
ggsave("pesudotime.group.png", plot4)
```

```{r}
#State轨迹分布图
plot1 <- plot_cell_trajectory(mycds, color_by = "State")
ggsave("State.pdf", plot = plot1, width = 6, height = 5)
ggsave("State.png", plot = plot1, width = 6, height = 5)
##Cluster轨迹分布图
plot2 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters")
ggsave("Cluster.pdf", plot = plot2, width = 6, height = 5)
ggsave("Cluster.png", plot = plot2, width = 6, height = 5)
##Pseudotime轨迹图
plot3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")
ggsave("Pseudotime.pdf", plot = plot3, width = 6, height = 5)
ggsave("Pseudotime.png", plot = plot3, width = 6, height = 5)
##合并作图
plotc <- plot1|plot2|plot3
ggsave("Combination.pdf", plot = plotc, width = 10, height = 3.5)
ggsave("Combination.png", plot = plotc, width = 10, height = 3.5)
```

```{r}
p1 <- plot2+facet_wrap( ~ orig.ident)
p2 <- plot2+facet_wrap( ~ seurat_clusters)
grid.arrange(p1, p2,nrow=2)
```


```{r}
diff_test_res <- differentialGeneTest(data1, fullModelFormulaStr = "~orig.ident")
# Select genes that are significant at an FDR < 10%
sig_genes <- subset(diff_test_res, qval < 0.1)
dim(sig_genes)
```


#6. 按组染色
#6.1 按组染色组合
```{r}
umap3 = data1@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% cbind(tx = data1@meta.data$orig.ident) %>% 
  mutate(group = substring(tx,1,2))

## 用 ggplot 绘图
p <- ggplot(umap3, aes(x = UMAP_1, y = UMAP_2, color = group), pt.size=1) + 
  geom_point(size = 1, alpha = 0.6) + 
  scale_color_manual(values = c("#E64B35FF","#25859B","#00A087B2"))+
  gait.errorcheck()+
  theme_void()

ggsave(file = "umap.by_group.pdf", width = 5, height = 5)
ggsave(file = "umap.by_group.png", width = 5, height = 5, units = "in", dpi = 300)
```


#6.2 按组分p
```{r}
#substring-截取字符
## 用 ggplot 绘图
p <- ggplot(umap3, aes(x = UMAP_1, y = UMAP_2, color = group), pt.size=1) + 
  geom_point(size = 1, alpha = 0.6) + 
  scale_color_manual(values = c("#E64B35FF","#25859B","#00A087B2"))+
  gait.errorcheck()+
  facet_wrap( ~ group)+
  theme_void()

ggsave(file = "umap.by_each_group.pdf", width = 10, height = 3.5)
ggsave(file = "umap.by_each_group.png", width = 10, height = 3.5, units = "in", dpi = 300)
```

#6.3 按样本染色 按组分p
```{r}
colourCount = length(unique(umap3$tx))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))


#组内一致性
p <- ggplot(umap3, aes(x = UMAP_1, y = UMAP_2, color = tx), pt.size=1) + 
  geom_point(size = 1, alpha = 0.6) + 
  scale_color_manual(values = getPalette(colourCount))+
  gait.errorcheck()+
  facet_wrap( ~ group)+
  theme_void()

ggsave(file = "umap.by_each_sample.pdf", width = 10, height = 3.5)
ggsave(file = "umap.by_each_sample.png", width = 10, height = 3.5, units = "in", dpi = 300)
```



#6. 寻找经典的marker基因，染色
```{r}
cellmarkers <- read.delim("filtered cell markers.txt", quote= "") $markers %>% as.character()
p <- FeaturePlot(data1, features = cellmarkers, pt.size = 0.5, reduction = "umap")
ggsave(file = "umap.cell.markers.png",p, width = 30, height = 30, limitsize = FALSE)
```

#6.1 以M1/M2 marker gene 染色
```{r}
genelist <- data1@assays$RNA@data@Dimnames[[1]]
geneset.IL <- genelist[grep("IL",genelist)]
for(i in seq(1, 126, by = 20 )){
  p <- FeaturePlot(data1, features = geneset.IL[i:(20+i)], pt.size = 0.3, reduction = "umap")
  ggsave(file = paste0("./cell.marker染色/IL.",i,".png"),p, width = 20, height = 20,units = "in", dpi = 300, limitsize = FALSE)
}

geneset.IL.filter <- c("IL6R","IL12RB2","SKIL","IL7R","IL6ST","CARMIL1","IL13RA1","IL15RA","IL18")
p <- FeaturePlot(data1, features = geneset.IL.filter, pt.size = 0.3, reduction = "umap")
ggsave(file = paste0("./cell.marker染色/IL.filter.png"),p, width = 20, height = 20,units = "in", dpi = 300, limitsize = FALSE)

geneset.TGF <- genelist[grep("TGF",genelist)]
p <- FeaturePlot(data1, features = geneset.TGF, pt.size = 0.3, reduction = "umap")
ggsave(file = paste0("./cell.marker染色/TGFr.png"),p, width = 20, height = 20,units = "in", dpi = 300, limitsize = FALSE)

geneset.TNF <- genelist[grep("TNF",genelist)]
p <- FeaturePlot(data1, features = geneset.TNF, pt.size = 0.3, reduction = "umap")
ggsave(file = paste0("./cell.marker染色/TNFr.png"),p, width = 20, height = 40,units = "in", dpi = 300, limitsize = FALSE)

geneset.IFN <- genelist[grep("IFN",genelist)]
p <- FeaturePlot(data1, features = geneset.IFN, pt.size = 0.3, reduction = "umap")
ggsave(file = paste0("./cell.marker染色/IFN.png"),p, width = 20, height = 40,units = "in", dpi = 300, limitsize = FALSE)
```


#7. cluster marker heatmap
```{r}
#加载cluster marker
memory.limit(1000000)
clus.marker <- read.delim("GraphClust.AllMarkerGenes.txt") %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_logFC)

pdf(file = "heatmap cluster.pdf", width = 9, height = 4)
png(file = "heatmap cluster.png", width = 15, height = 8, units = "in", res = 300)
DoHeatmap(data1, features = clus.marker$allmarkers.gene)+ NoLegend()
dev.off()
```
#8. bubbleplot cluster
```{r}
#bubble plot，由上往下看，顺序要反过来
# my_levels <- read.delim("culs level bubble.txt", quote= "") $order %>% as.character()
# Idents(data) <- factor(Idents(data), levels= my_levels)
DotPlot(data1, features = clus.marker$allmarkers.gene, cols = c("#FEE0D2", "#CB181D"), scale.min = 10) + 
  scale_x_discrete("")+
  scale_y_discrete("")+
  RotatedAxis()
ggsave(file = "Bub cluster.png", width = 10, height = 8, limitsize = FALSE)
ggsave(file = "Bub cluster.pdf", width = 10, height = 8, limitsize = FALSE)
```

#9. S3M染色
#9.1 cluster染色(用不上)
```{r}
umap1 = data1@reductions$umap@cell.embeddings %>%
      as.data.frame() %>% cbind(tx = data1@meta.data$seurat_clusters)

#颜色梯度
Palette.S3M <- c("grey", "grey", "grey", "#C4625D", "grey", "grey", "grey", "grey")

## 用 ggplot 绘图
ggplot(umap1, aes(x = UMAP_1, y = UMAP_2, color = tx), pt.size=1) + 
  geom_point(size = 1, alpha = 0.6) + 
  scale_color_manual(values = Palette.S3M)+
  theme_void()
ggsave(file = "./S3M/umap.by_cluster.pdf", width = 4, height = 4)
ggsave(file = "./S3M/umap.by_cluster.png", width = 4, height = 4, units = "in", dpi = 300)
```

#9.2 只S3M染色，按组分p
```{r}
umap3 = umap1 %>% cbind(sample = data1@meta.data$orig.ident) %>% 
  mutate(group = substring(sample,1,2))

## 用 ggplot 绘图
p1 <- ggplot(umap3%>%filter(group == "De"), aes(x = UMAP_1, y = UMAP_2, color = tx), pt.size=1) + 
  geom_point(size = 1, alpha = 0.6) + 
  scale_color_manual(values = c("grey", "grey", "grey", "#E64B35FF", "grey", "grey", "grey", "grey"))+
  theme_void()

p2 <- ggplot(umap3%>%filter(group == "HC"), aes(x = UMAP_1, y = UMAP_2, color = tx), pt.size=1) + 
  geom_point(size = 1, alpha = 0.6) + 
  scale_color_manual(values = c("grey", "grey", "grey", "#25859B", "grey", "grey", "grey", "grey"))+
  theme_void()

p3 <- ggplot(umap3%>%filter(group == "Re"), aes(x = UMAP_1, y = UMAP_2, color = tx), pt.size=1) + 
  geom_point(size = 1, alpha = 0.6) + 
  scale_color_manual(values = c("grey", "grey", "grey", "#00A087B2", "grey", "grey", "grey", "grey"))+
  theme_void()
p1|p2|p3

ggsave(file = "./S3M/umap.by_each_group.pdf", width = 10, height = 3.5)
ggsave(file = "./S3M/umap.by_each_group.png", width = 10, height = 3.5, units = "in", dpi = 300)
```

#9.3 
```{r}
my_comparisons <- list(c("De", "HC"), c("De", "Re"), c("HC", "Re"))
p = ggboxplot(cell.info %>% filter(seurat_clusters == 3) , x="group", y="cell.per", color = "group",
               palette = c("#FC4E07","#00AFBB", "#E7B800"), 
               add = "jitter")+
    theme(legend.title=element_blank())+
  stat_compare_means(method = "wilcox.test", comparisons = my_comparisons)+
  stat_compare_means(label.y = 1)
```


#10 S03M vs. 小胶质DEG overlap基因染色和小提琴
```{r} 
data1@meta.data <- data1@meta.data %>% mutate(group = substr(orig.ident, 1,2))
genelist.56 <- read.delim("S3M/S03M vs MICdeg.txt")$DEGs
genelist.15 <- read.delim("S3M/S03M vs MICdeg filtered.txt")$DEGs

VlnPlot(data1, features = genelist.56[41:56], group.by = "group", pt.size =0)
VlnPlot(data1, features = genelist.15, cols = c("#E64B35E5", "#4DBBD5E5", "#00A087E5"), ncol = 6, group.by = "group", pt.size = 0)

ggsave(file = "./S3M/vln.56gene.filtered.pdf", width = 12, height = 10)

Idents(data1) <- "group"
micro_deg_res <- FindMarkers(data1, ident.1 = "De", ident.2= "Re")
micro_deg_HC <- FindMarkers(data1, ident.1 = "De", ident.2= "HC")
write.xlsx(micro_deg_HC,"microglia.dep-hc.xlsx",row.names = T)
write.xlsx(micro_deg_res,"microglia.dep-res.xlsx",row.names = T)
```


#10.2 S03M vs. 小胶质DEG overlap基因 精神疾病注释
```{r}
library(psygenet2r)
library(disgenet2r)
disgenet_api_key <- get_disgenet_api_key(email = "470355362@qq.com",password = "wujing301" )

psy.anno.56 <- psygenetGene(genelist.56, database = "ALL", verbose = FALSE)

pdf("56gene.output.pdf")
plot(psy.anno.56)
plot(psy.anno.56, 'GDCA network')
plot(psy.anno.56, 'heatmapGenes')
geneAttrPlot(psy.anno.56, type = "evidence index")
```










