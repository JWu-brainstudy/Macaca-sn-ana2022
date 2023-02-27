#加载包
library(Seurat)
library(dplyr)
library(tidyverse)
#加载数据
snRNA1 <- readRDS("macacards/snRNA.1.rds")
#降维
snRNA1 <- ScaleData(snRNA1, features = (rownames(snRNA1)))
snRNA1 <- RunPCA(snRNA1, features = VariableFeatures(snRNA1),seed.use=3)
DimPlot(snRNA1, reduction = "pca", group.by="stim")
snRNA1 <- JackStraw(snRNA1,reduction = "pca", dims=50)
snRNA1 <- ScoreJackStraw(snRNA1,dims = 1:50)
ElbowPlot(snRNA1, ndims=20, reduction="pca")
JackStrawPlot(snRNA1, dims = 1:20)
#聚类
data1 <- FindNeighbors(snRNA1, dims = 1:30) %>% 
	RunUMAP(dims = 1:10) %>% 
	FindClusters(resolution = 0.8)
#save
saveRDS(data1,file = "snRNA.2.rds")