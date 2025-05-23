library(ggplot2)
library(ggnewscale)
library(cowplot)
library(dplyr)
library(ggpubr)
library(scCustomize)

d<-read.csv("/Users/amistry/Library/CloudStorage/OneDrive-uoflhealth/Atlas/Atlas_data.csv", header=T, row.names = 1, stringsAsFactors = T)
d<-d[d$Diagnosis!="#N/A",]
d$Diagnosis<-factor(d$Diagnosis, levels = c("Papilloma", "Angiocentric", "Diffuse Glioma-WT", "G34 Mutant", "IDH Mutant", "IDH Mut-1p19q del", "Midline",
                                    "HGNET-MN1", "PA", "PXA", "SEGA", "AT/RT-MYC", "AT/RT-SHH", "AT/RT-TYR", "ETMR", "HGNET-BCOR","MB-GP3",
                                    "MB-GP4","MB-SHH","MB-WNT", "Neuroblastoma", "NB-FOXR2", "Retinoblastoma", "MPE", "PFA", "PFB", "RELA", "SE", "YAP", "Germ Cell Tumor",
                                    "DIG", "Ganglioglioma", "PCNSL", "Uveal Melanoma", "EFT-CIC", "Meningioma", "SFT", "Ganglioneuroma",
                                    "MPNST", "Neurofibroma", "PH/PG", "Ada. Cranio.", "PitNET", "Breast", "Lung",
                                    "Cerebellum", "Cerebellum, Fetal", "Choroid Plexus", "CNS", "Nerve Ganglia", "Peri. Nerve",
                                    "Pituitary", "Retina", "Retina, Fetal"))

colors = c(DiscretePalette_scCustomize(num_colors = 45, palette = "varibow", shuffle_pal = F),"black","grey10","grey20","gray80","grey40","grey50","grey60","grey70","grey30")

colors<-c("#CC985C",  
          "#994550", "#FFD426", "#CC0036", "#FF00CC","#FF267D", "#99752E", 
          "#939945", "#CCBE00", "#FF7373", "#92CC1F",
          "#529900", "#944DFF", "#63CC3D", "#289917", "#992E67", "#5CCC6B", "#2E994A", "#26FF7D", "#00CC6D", "#DBFF4D", "#1FCCC0", "#4DFFDB", 
          "#008F99", "#73E3FF", "#3D9CCC", "#175C99", "#0066FF", "#C23DCC",
          "#2E3C99", 
          "#2626FF", "#1B00CC", 
          "#5B4599", 
          "#00FF00", 
          "#7B1FCC", "#45997D", "#E373FF", 
          "#5C7ACC", "#991790", "#CC503D", "#ABFF73",
          "#CC5CA7", "#660099", 
          "#993A17", "#FF6600", "black", "grey10", "grey20", "gray80", "grey40", "grey50", "grey60", "grey70", "grey30")

{
  p<-ggplot()+geom_point(data=d, aes(tsne1,tsne2,color=factor(Diagnosis)), size=0.75)+theme_pubr()+
    #geom_text(data = lab, aes(x=x,y=y,label=Diagnosis), size=3)+
    annotate("text", x=-94, y=80, label= "t-SNE Dimensionality Reduction\nBefore SVA Adjustment\n(5264 samples)", hjust = 0, fontface = "bold") + 
    #annotate("text", x=-94, y=-75, size=3, label= "Angio. = Angiocentric\nBrMet = Breast Cancer Metastasis\nGCT = Germ Cell Tumor\nGN = Ganglioneuroma\nSEGA = Subependymal Giant Cell Astrocytoma\nUM = Uveal Melanoma", hjust = 0) + 
    scale_color_manual(values=colors)+
    theme(legend.position = "none")+rremove("axis")+rremove("axis.text")+rremove("ticks")+rremove("xylab")+
    theme(panel.border = element_rect(colour = "black", fill=NA))
  
  p1<-ggplot()+
    geom_point(data=d[d$Diagnosis.class=="Choroid Plexus",],aes(tsne1,tsne2,color=Diagnosis),size=1)+
    scale_color_manual(values=colors[1], name="Choroid Plexus")+
    guides(color = guide_legend(override.aes=list(shape = 15, size=5), order = 1))+
    new_scale_color() +
    
    geom_point(data=d[d$Diagnosis.class=="Diffuse Glioma",],aes(tsne1,tsne2,color=Diagnosis),size=1)+
    scale_color_manual(values=colors[2:7], name="Diffuse Glioma")+
    guides(color = guide_legend(override.aes=list(shape = 15, size=5, ncol=2), order = 2))+
    new_scale_color() +
    
    geom_point(data=d[d$Diagnosis.class=="Discrete Glioma",],aes(tsne1,tsne2,color=Diagnosis),size=1)+
    scale_color_manual(values=colors[8:11], name="Discrete Glioma")+
    guides(color = guide_legend(override.aes=list(shape = 15, size=5), order = 3))+
    new_scale_color() +
    
    geom_point(data=d[d$Diagnosis.class=="Embryonal",],aes(tsne1,tsne2,color=Diagnosis),size=1)+
    scale_color_manual(values=colors[12:23], name="Embryonal")+
    guides(color = guide_legend(override.aes=list(shape = 15, size=5), order = 4))+
    
    theme_classic()+theme(legend.title=element_text(face="bold"), legend.justification = "top")
  
  p1l <- get_legend(p1)
  
  p2<-ggplot()+
    geom_point(data=d[d$Diagnosis.class=="Ependymal",],aes(tsne1,tsne2,color=Diagnosis),size=1)+
    scale_color_manual(values=colors[24:29], name="Ependymal")+
    guides(color = guide_legend(override.aes=list(shape = 15, size=5), order = 1))+
    new_scale_color() +
    
    geom_point(data=d[d$Diagnosis.class=="Germ Cell",],aes(tsne1,tsne2,color=Diagnosis),size=1)+
    scale_color_manual(values=colors[30], name="Germ Cell")+
    guides(color = guide_legend(override.aes=list(shape = 15, size=5), order = 2))+
    new_scale_color() +
    
    geom_point(data=d[d$Diagnosis.class=="Glioneuronal",],aes(tsne1,tsne2,color=Diagnosis),size=1)+
    scale_color_manual(values=colors[31:32], name="Glioneuronal")+
    guides(color = guide_legend(override.aes=list(shape = 15, size=5), order = 3))+
    new_scale_color() +
    
    geom_point(data=d[d$Diagnosis.class=="Hematolymphoid",],aes(tsne1,tsne2,color=Diagnosis),size=1)+
    scale_color_manual(values=colors[33], name="Hematolymphoid")+
    guides(color = guide_legend(override.aes=list(shape = 15, size=5), order = 4))+
    new_scale_color() +
    
    geom_point(data=d[d$Diagnosis.class=="Melanocytic",],aes(tsne1,tsne2,color=Diagnosis),size=1)+
    scale_color_manual(values=colors[34], name="Melanocytic")+
    guides(color = guide_legend(override.aes=list(shape = 15, size=5), order = 5))+
    new_scale_color() +
    
    geom_point(data=d[d$Diagnosis.class=="Mesenchymal",],aes(tsne1,tsne2,color=Diagnosis),size=1)+
    scale_color_manual(values=colors[35:37], name="Mesenchymal")+
    guides(color = guide_legend(override.aes=list(shape = 15, size=5), order = 6))+
    
    theme_classic()+theme(legend.title=element_text(face="bold"),legend.justification = "top")
  
  p2l <- get_legend(p2)
  
  p3<-ggplot()+
    geom_point(data=d[d$Diagnosis.class=="Nerve",],aes(tsne1,tsne2,color=Diagnosis),size=1)+
    scale_color_manual(values=colors[38:41], name="Nerve")+
    guides(color = guide_legend(override.aes=list(shape = 15, size=5), order = 1))+
    new_scale_color() +
    
    geom_point(data=d[d$Diagnosis.class=="Sellar",],aes(tsne1,tsne2,color=Diagnosis),size=1)+
    scale_color_manual(values=colors[42:43], name="Sellar")+
    guides(color = guide_legend(override.aes=list(shape = 15, size=5), order = 2))+
    new_scale_color() +
    
    geom_point(data=d[d$Diagnosis.class=="Metastasis",],aes(tsne1,tsne2,color=Diagnosis),size=1)+
    scale_color_manual(values=colors[44:48], name="Metastasis")+
    guides(color = guide_legend(override.aes=list(shape = 15, size=5), order = 3))+
    new_scale_color() +
    
    geom_point(data=d[d$Diagnosis.class=="Non-tumor",],aes(tsne1,tsne2,color=Diagnosis),size=1)+
    scale_color_manual(values=c("black","grey10","grey20","gray80","grey40","grey50","grey60","grey70","grey30"), name="Non-tumor")+
    guides(color = guide_legend(override.aes=list(shape = 15, size=5), order = 4))+
    new_scale_color() +
    
    theme_classic()+theme(legend.title=element_text(face="bold"),legend.justification = "top")
  
  p3l <- get_legend(p3)
  
  cowplot::plot_grid(p,p1l,p2l,p3l, ncol=4, rel_widths = c(6,1,1,1), align = 'hv')
}
#1200 x 717 svg

