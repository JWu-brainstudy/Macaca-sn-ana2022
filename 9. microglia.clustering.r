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
library(hrbrthemes)
library(cowplot)
library(xlsx)
library(reshape)
library(rstatix)
library(ggsignif)

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
write.csv(cell.info,"小胶相对丰度.csv")


cell<-read.csv("小胶相对丰度.csv") %>% mutate(seurat_clusters = paste("Mic",seurat_clusters))

# 分组计算P值
stat.test <- cell %>% group_by(seurat_clusters)%>%
        wilcox_test(cell.per~group)%>% add_significance()

p.m3 <- ggplot(cell %>% filter(seurat_clusters == "Mic 3"),
               aes(group,cell.per,fill=group))+
        geom_boxplot(linetype = "dashed", outlier.alpha = 0, alpha=0)+ #全体虚线
        stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.alpha = 0, alpha=0.5)+ #实线盖住中间框
        stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.4)+ #上部实线
        stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width=0.4)+ #上部实线
        theme_classic(base_family="Arial")+ylim(0,0.4) +
        scale_fill_manual(values = c("HC"="#46B1C9", "De"="#DA4C36", "Re"="#229982"))+
        geom_dotplot(binaxis='y', stackdir='center', stackratio=1, dotsize=0.5) + 
        scale_x_discrete(limits=c("De", "Re", "HC"))

ggsave("M03主图.pdf", p.m3, device = "pdf", width = 5,height = 10)

p.all.2023 <- ggplot(cell, aes(group,cell.per,fill=group))+
        geom_boxplot(linetype = "dashed", outlier.alpha = 0, alpha=0)+ #全体虚线
        stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..), outlier.alpha = 0, alpha=0.1)+ #实线盖住中间框
        stat_boxplot(geom = "errorbar", aes(ymin = ..ymax..), width=0.4)+ #上部实线
        stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..), width=0.4)+ #上部实线
        theme_ipsum(base_family="Arial")+ylim(0,0.4) +
        scale_fill_manual(values = c("HC"="#46B1C9", "De"="#DA4C36", "Re"="#229982"))+
        geom_signif(comparisons = list(c("De","HC"),c("De","Re"),c("HC","Re")),
              map_signif_level=c("***"=0.001,"**"=0.01,"*"=0.05),
              textsize=3,test=wilcox.test,step_increase=0.25)+
        geom_dotplot(binaxis='y', stackdir='center', stackratio=1.1, dotsize=0.5) +
        scale_x_discrete(limits=c("De", "Re", "HC"))+
        facet_wrap(~seurat_clusters ,scales="free", ncol = 8)

ggsave("subcluster.2023.2.pdf", p.all.2023, device = "pdf", width = 16,height = 4)


#小胶质细胞相关marker
color <- c("#EAEAEA", "#e41a1c", "#3f007d")
p_iq <- FeaturePlot(data1, features = c("IQGAP2", "FYN", "PDE7A",  "ARHGEF3","BCL11B"), cols = color, pt.size = 0.5, reduction = "umap", ncol = 5)+ theme_void()
ggsave(file = "markers.IQGAP2.png", width = 20, height = 4, limitsize = FALSE)



# 按组染色
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


# 按组分p
p <- ggplot(umap3, aes(x = UMAP_1, y = UMAP_2, color = group), pt.size=1) + 
  geom_point(size = 1, alpha = 0.6) + 
  scale_color_manual(values = c("#E64B35FF","#25859B","#00A087B2"))+
  gait.errorcheck()+
  facet_wrap( ~ group)+
  theme_void()
ggsave(file = "umap.by_each_group.pdf", width = 10, height = 3.5)
ggsave(file = "umap.by_each_group.png", width = 10, height = 3.5, units = "in", dpi = 300)


# Mic03
```{r}
umap3 = umap1 %>% cbind(sample = data1@meta.data$orig.ident) %>% 
  mutate(group = substring(sample,1,2))

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

# boxplot 
my_comparisons <- list(c("De", "HC"), c("De", "Re"), c("HC", "Re"))
p = ggboxplot(cell.info %>% filter(seurat_clusters == 3) , x="group", y="cell.per", color = "group",
               palette = c("#FC4E07","#00AFBB", "#E7B800"), 
               add = "jitter")+
    theme(legend.title=element_blank())+
  stat_compare_means(method = "wilcox.test", comparisons = my_comparisons)+
  stat_compare_means(label.y = 1)

# S03M vs. 小胶质DEG overlap基因染色和小提琴
data1@meta.data <- data1@meta.data %>% mutate(group = substr(orig.ident, 1,2))
genelist.15 <- read.delim("S3M/S03M vs MICdeg filtered.txt")$DEGs
VlnPlot(data1, features = genelist.15, cols = c("#E64B35E5", "#4DBBD5E5", "#00A087E5"), ncol = 6, group.by = "group", pt.size = 0)

ggsave(file = "./S3M/vln.56gene.filtered.pdf", width = 12, height = 10)


# S03M vs. 小胶质DEG overlap基因 精神疾病注释
library(psygenet2r)
library(disgenet2r)
disgenet_api_key <- get_disgenet_api_key(email = "470355362@qq.com",password = "******" )

psy.anno.56 <- psygenetGene(genelist.56, database = "ALL", verbose = FALSE)

pdf("56gene.output.pdf")
plot(psy.anno.56)
plot(psy.anno.56, 'GDCA network')
plot(psy.anno.56, 'heatmapGenes')
geneAttrPlot(psy.anno.56, type = "evidence index")
