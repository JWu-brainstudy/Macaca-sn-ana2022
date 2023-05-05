library(devtools)
# install_bitbucket("ibi_group/disgenet2r")
library(disgenet2r)
library(psygenet2r)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

disgenet_api_key <- get_disgenet_api_key(email = "470355362@qq.com",password = "*******" )

degs <- read.delim("DEGS.txt",sep = "\t")
query_degs <- as.character(degs$Gene) %>% unique(fromLast = TRUE)
data1 <- gene2disease( gene = query_degs,score =c(0.2, 1), verbose = TRUE, api_key = disgenet_api_key)
plot(data1,class="Heatmap",limit=100)
ggsave("deg2disease.heatmap.png", width = 12, height = 12)
ggsave("deg2disease.heatmap.pdf", width = 12, height = 12)
plot(data1,class="DiseaseClass")
ggsave("deg2disClass.pdf", width = 12, height = 12)
ggsave("deg2disClass.png", width = 12, height = 12)

#Gene disease Heatmap

dis.erich <- read.delim("deg2disease.filtered.txt") %>% filter(genes > 7)

p.orign <- ggplot(dis.erich, aes(gene_symbol, disease_name)) + 
  geom_tile(aes(fill = score))+
  scale_fill_gradient2(low = "#FFF5F0",high = "#A50F15") +
  theme_classic() + 
  theme(axis.title.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.title.y=element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 12)) + 
  #labs(fill =paste0(" * p < 0","\n\n","** p > 0","\n\n","Correlation")) +  
  scale_x_discrete(position = "top") 
print(p.orign)
                   
ggsave("dis2gene.png", width = 10, height = 10)
ggsave("dis2gene.pdf", width = 10, height = 10)

data2 <- psygenetGeneSentences(query_degs, database = "ALL", verbose = FALSE)

pdf("output.pdf")

plot(data2)
plot(data2, "GDCA network")
plot(data2, "heatmapGenes")
geneAttrPlot(data2, type = "evidence index")
function.disease <- pantherGraphic(query_degs, "ALL")
pantherGraphic(query_degs[1:5], "ALL")


