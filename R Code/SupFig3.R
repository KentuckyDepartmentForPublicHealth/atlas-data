library(ggplot2)
library(ggpubr)

d<-read.csv("/Users/amistry/Library/CloudStorage/OneDrive-uoflhealth/Atlas/tsne_74_5000.csv", header=T, row.names = 1, stringsAsFactors = T)

ggplot()+geom_point(data=d, aes(V1,V2,color=Sample), size=0.75)+theme_pubr()+
  annotate("text", x=-94, y=87, label= "t-SNE Dimensionality Reduction\n(7375 samples)", hjust = 0, fontface = "bold") + 
  scale_color_manual(values=c("grey","black","red"))+
  theme(legend.position = "right")+rremove("axis")+rremove("axis.text")+rremove("ticks")+rremove("xylab")+
  theme(panel.border = element_rect(colour = "black", fill=NA))
#975x670