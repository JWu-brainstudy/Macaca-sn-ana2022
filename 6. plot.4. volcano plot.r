
library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggrepel)
library(RColorBrewer)
library(xlsx)
library(dplyr)
library(ggsci)
library(scales)


#1加载数据
data <-read.xlsx("逆转基因作图.xlsx",1)%>%
  mutate(color = paste(pattern,cell,sep="_"))

#2 挑选颜色
show_col(pal_lancet(alpha =1)(8))
show_col(pal_lancet(alpha =0.6)(8))
pal_lancet(alpha =1)(8)
pal_lancet(alpha =0.6)(8)
show_col(pal_nejm(alpha =1)(8))
show_col(pal_nejm(alpha =0.6)(8))
pal_nejm(alpha =1)(8)
pal_nejm(alpha =0.6)(8)

volcano <- ggplot(data, aes(avg_logFC, -1*log10(p_val_adj)))+ 
  geom_point(aes(color = color), size = 2, alpha = 0.4)+
  geom_label_repel(aes(label = Genename, color = color), 
                   size = 3, max.overlaps = 20,fontface = 'italic') + 
  labs(x = expression(log10(FC)), y = expression(-log[10](FDR))) +
  scale_color_manual(values = c("#ED0000FF","#00468BFF","#FFDC91FF","#42B540FF","#ED000099","#0099B499","#925E9F99","#00468B99","#FFDC9199","#42B54099")) + 
  geom_vline(xintercept = 0, linetype = 4, size = 1, color = "grey")+
  ggtitle("DEG")+
  theme_classic()+
  xlim(-1.2, 1.2)+
  theme(axis.line = element_line(size = 1, color = "black"))+
  theme(axis.text = element_text(size=10, color='#003087'))

print(volcano)
ggsave("主图逆转基因.pdf",volcano, width = 10, height = 8)
ggsave("./主图逆转基因.png",volcano, width = 10, height = 8)
