library(dplyr)
library(openxlsx)
library(Seurat)

data <- readRDS("macacards/snRNA.2.rds")
Idents(data) = data$Celltype
microglia <- SplitObject(object, split.by = "ident")[["Mic"]]
saveRDS(microglia, "microglia.subset.rds")

