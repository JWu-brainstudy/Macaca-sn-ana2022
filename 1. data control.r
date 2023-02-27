#加载包
library(Seurat)
library(dplyr)
library(tidyverse)
#导入数据
samples=list.files("GSE201687/")
samples
dir <- file.path('./GSE201687',samples)
names(dir) <- samples
#构建rds对象
counts <- Read10X(data.dir = dir)
snRNA1 = CreateSeuratObject(counts, min.cells = 3, min.features = 100)
#数据质控
snRNA1[["percent.mt"]] <- PercentageFeatureSet(snRNA1, pattern = "^mt-")
pctMT = 5
snRNA1 <- subset(snRNA1, subset = percent.mt < pctMT)
VlnPlot(snRNA1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0.05)
#可变基因标记
snRNA1 <- FindVariableFeatures(snRNA1, selection.method = "vst", nfeatures = 2000) 
top10 <- head(VariableFeatures(snRNA1), 10) 
top10
plot1 <- VariableFeaturePlot(snRNA1)
#save
saveRDS(snRNA1,file = "macacards/snRNA.1.rds")

