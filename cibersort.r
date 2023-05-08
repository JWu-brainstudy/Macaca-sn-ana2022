#Set the path
setwd("../macacards/")

#Install Package
library(devtools)
devtools::install_github("Moonerss/CIBERSORT")

#Load Packages
library(dplyr)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggrepel)
library(CIBERSORT)

#Activate function
source("CIBERSORT.R")

#Load microglia subclustered
scRNA <- readRDS("microglia.subset.rds") 

#Cluster marker genes screening
X <- FindAllMarkers(scRNA) %>% 
  filter(p_val_adj < 0.05) %>%
  filter(avg_log2FC > 0.25)

#Make a cluster coefficient table
text = scRNA[['RNA']]@scale.data
text = text[rownames(text) %in% !duplicated(X$gene),]
write.table(text,"text.txt",sep = "\t",col.names = T,row.names = T) 

sig <- AverageExpression(scRNA)[[1]]  
write.table(sig,"sig.txt",sep = "\t",col.names = T,row.names = T) 

#Cibersort Analyze
res_cibersort=CIBERSORT("sig.txt","test.txt", perm=1000, QN=TRUE) 
save(res_cibersort,file = "res_cibersort.Rdata")  

#Cibersort data processing
rm(list=ls())
load("res_cibersort.Rdata")
res_cibersort <- res_cibersort[,1:22] #The first 22 columns were taken as cell abundance data
ciber.res <- res_cibersort[,colSums(res_cibersort) > 0] #Remove cells with all 0 abundance

#Barplot plot
mycol <- ggplot2::alpha(rainbow(ncol(ciber.res)), 0.7) 
par(bty="o", mgp = c(2.5,0.3,0), mar = c(2.1,4.1,2.1,10.1),tcl=-.25,las = 1,xpd = F)
a = barplot(as.matrix(t(ciber.res)),
            border = NA, 
            names.arg = rep("ST region",nrow(ciber.res)),
            yaxt = "n", 
            ylab = "Propotion", 
            col = mycol) 
axis(side = 2, at = c(0,0.2,0.4,0.6,0.8,1), 
     labels = c("0","0.2","0.4","0.6","0.8","1"))
legend(par("usr")[2]-20, 
       par("usr")[4], 
       legend = colnames(ciber.res), 
       xpd = T,
       fill = mycol,
       cex = 0.8, 
       border = NA, 
       y.intersp = 1,
       x.intersp = 0.2,
       bty = "n")
dev.off()
