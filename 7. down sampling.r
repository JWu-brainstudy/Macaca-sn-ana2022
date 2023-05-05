library(dplyr)
library(openxlsx)
library(Seurat)

cell.list <- read.delim("summarycell_celltype.txt")

summarise(group_by(cell.list, cellType), n())

#查看别人的抽样比例，自闭症选择1900细胞, 是因为IN-SV2C有1914个细胞

ASD.cell.list <- read.xlsx("D:/singlecellpapers/autism/aav8130_Data-S2.xlsx")[1:4]
summarise(group_by(ASD.cell.list, cluster), n())

#去掉undentified后，按3700细胞随机取样,重复10次
downsampling.list <- list(c())
for(i in 1:10){
  downsampling.list[[i]] <- filter(cell.list, CellType != "Undefined") %>%
    group_by(CellType) %>%
    sample_n(3700)
}

#测试抽样结果
intersect(downsampling.list[[2]],downsampling.list[[4]])
summarise(group_by(downsampling.list[[2]], CellType), n())

#去掉undentified和Fibro&SMC&Endo后，按3700细胞随机取样,重复10次
downsampling.list <- list(c())
for(i in 1:10){
  downsampling.list[[i]] <- filter(cell.list, CellType %in% c("Astrocyte", "Excitatory neuron", "Interneuron", "Microglia", "Oligodendrocyte", "OPC")) %>%
    group_by(CellType) %>%
    sample_n(3700)
}

filename <- paste("celllist.rep",1:10,".txt", sep = "")
for(i in 1:10){
  write.table(downsampling.list[[i]], filename[i], row.names = F, quote = F, sep = "\t")
}

#repeat 10 times,celllist.rep1.txt to celllist.rep10.txt
data <- readRDS("macacards/snRNA.2.rds")
Idents(data) = data$Group
cell.list <- read.delim2("celllist.rep1.txt", sep = "\t", quote = F)
data <- subset(data, cells = cell.list$Cell)
DEG.1 <- FindMarkers(data, ident.1 = Dep, ident.2 = HC, group.by = Celltype)
n.deg <- DEG.1 %>% group_by(Celltype) %>% summarise(n.deg = n())

