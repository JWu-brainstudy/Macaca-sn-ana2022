#Set the path
setwd("../macacards/")

#Load Packages
library(Seurat)
library(openxlsx)
library(dplyr)
library(ggsci)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(rstatix)
library(reshape)
library(WGCNA)

#Load all celltypes datasets 
#Extracting Int expression matrix
data1 <- readRDS("snRNA.2.rds")
Idents(data1) = data1$Celltype
x1 <- SplitObject(data1,split.by = "ident")[["Int"]]
x <- AverageExpression(x1,group.by = "orig.ident")$RNA 


#Load phenotypic data
trait <- read.csv("../traut.csv",header = T)[,5:11]

#Missing value handling
gsg = goodSamplesGenes(x, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
   if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(x)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(x)[!gsg$goodSamples], collapse = ", ")));
   x = x[gsg$goodSamples, gsg$goodGenes]
}

#Sample clustering
sampleTree = hclust(dist(x), method = "average");
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
ggsave("sampleClustering_Int.png",height = 9,width = 12)

#Soft threshold calculation
powers = c(c(1:20))
sft = pickSoftThreshold(x, powerVector = powers, verbose = 5,networkType = "signed")
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.8;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
ggsave("soft-thresholding power_Int.png",height = 6,width = 12)

#Building network in one step
net = blockwiseModules(x, power = sft$powerEstimate,
                       TOMType = "unsigned", minModuleSize =300,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       networkType = "signed",
                       verbose = 3)
sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
ggsave("cluster Dendrogram_Int.png",height = 9,width = 12)

unmergedColors = labels2colors(net$unmergedColors)
mergedColors   = labels2colors(net$colors)
plotDendroAndColors(
  net$dendrograms[[1]],
  cbind(unmergedColors[net$blockGenes[[1]]], mergedColors[net$blockGenes[[1]]]),
  c("Dynamic Tree Cut" , "Merged colors"),
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)
ggsave("cluster-dynamic Dendrogram_Int.png",height = 9,width = 12)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
meta2module <- as.data.frame(moduleLabels) %>% cbind(mergedColors)
write.xlsx(meta2module, "meta2module_Int.xlsx", overwrite = T, rowNames = T)
nGenes = ncol(x);
nSamples = nrow(x);
MEs0 = moduleEigengenes(x, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, trait, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

#graphics
MDD.wgcna <- melt(moduleTraitCor) %>% 
  cbind(melt(moduleTraitPvalue))%>%
  .[,c(1:3,6)]
colnames(MDD.wgcna)=c("v1","v2","cor","P")
MDD.wgcna <- MDD.wgcna %>%
  mutate(text = case_when( 
    P >= 0.05 ~ paste(" "), 
    P < 0.05 & P >= 0.01 ~ paste(round(P, 3), " *"),
    P < 0.01 & P >= 0.001~ paste(round(P, 3), " **"),
    P < 0.001 ~ paste(round(P, 3), " ***")))
p.orign <- ggplot(MDD.wgcna, aes(v1, v2)) + 
  geom_tile(aes(fill = cor), colour = "grey", size = 0.5)+
  scale_fill_gradient2(low = "#377eb8",mid = "white",high = "#e41a1c") +
  geom_text(aes(label=text),col ="black",size = 1.5) +
  theme_minimal() + 
  theme(axis.title.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.title.y=element_blank(), 
        axis.text.x = element_text(angle = 45, size = 8), 
        axis.text.y = element_text(size = 8)) + 
  scale_x_discrete(position = "top") 
p.orign
ggsave("wgcna_Int.pdf", width = 12, height = 5)
ggsave("wgcna_Int.png", width = 12, height = 5)
