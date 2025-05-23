library(ggplot2)
library(ggnewscale)
library(cowplot)
library(dplyr)
library(ggpubr)
library(scCustomize)

d<-read.csv("/Users/amistry/Library/CloudStorage/OneDrive-uoflhealth/Atlas/Atlas_data.csv", header=T, row.names = 1, stringsAsFactors = T)
d$Diagnosis<-d$Final.diagnosis_bc
d$tsne1<-d$tsne1_bc
d$tsne2<-d$tsne2_bc
d<-d[d$Diagnosis!="-",]
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

lab<- d %>% group_by(Diagnosis) %>% summarize(x = median(tsne1), y = median(tsne2))
{
  lab$Diagnosis<-as.character(lab$Diagnosis)
  lab[which(lab$Diagnosis=="Germ Cell Tumor"),"Diagnosis"]<-"GCT"
  lab[which(lab$Diagnosis=="Lung"),"Diagnosis"]<-"Lung Met"
  lab[which(lab$Diagnosis=="Breast"),"Diagnosis"]<-"Breast Met"
  lab[which(lab$Diagnosis=="Angiocentric"),"Diagnosis"]<-"Angio."
  lab[which(lab$Diagnosis=="Ganglioglioma"),"Diagnosis"]<-"GG"
  lab$Diagnosis<-factor(lab$Diagnosis)
  lab[which(lab$Diagnosis=="Papilloma"),"x"]<-lab[which(lab$Diagnosis=="Papilloma"),"x"]+7.5
  lab[which(lab$Diagnosis=="Papilloma"),"y"]<-lab[which(lab$Diagnosis=="Papilloma"),"y"]+2.5
  lab[which(lab$Diagnosis=="Choroid Plexus"),"x"]<-lab[which(lab$Diagnosis=="Choroid Plexus"),"x"]+11
  lab[which(lab$Diagnosis=="Choroid Plexus"),"y"]<-lab[which(lab$Diagnosis=="Choroid Plexus"),"y"]-2.5
  lab[which(lab$Diagnosis=="YAP"),"x"]<-lab[which(lab$Diagnosis=="YAP"),"x"]+4
  lab[which(lab$Diagnosis=="YAP"),"y"]<-lab[which(lab$Diagnosis=="YAP"),"y"]+2.5
  lab[which(lab$Diagnosis=="PFB"),"x"]<-lab[which(lab$Diagnosis=="PFB"),"x"]-6
  lab[which(lab$Diagnosis=="PFB"),"y"]<-lab[which(lab$Diagnosis=="PFB"),"y"]+3
  lab[which(lab$Diagnosis=="SE"),"x"]<-lab[which(lab$Diagnosis=="SE"),"x"]-5.5
  lab[which(lab$Diagnosis=="SE"),"y"]<-lab[which(lab$Diagnosis=="SE"),"y"]+2
  lab[which(lab$Diagnosis=="MPE"),"x"]<-lab[which(lab$Diagnosis=="MPE"),"x"]-5.5
  lab[which(lab$Diagnosis=="RELA"),"x"]<-lab[which(lab$Diagnosis=="RELA"),"x"]+6
  lab[which(lab$Diagnosis=="RELA"),"y"]<-lab[which(lab$Diagnosis=="RELA"),"y"]+5
  lab[which(lab$Diagnosis=="Neuroblastoma"),"y"]<-lab[which(lab$Diagnosis=="Neuroblastoma"),"y"]+10
  lab[which(lab$Diagnosis=="Ganglioneuroma"),"y"]<-lab[which(lab$Diagnosis=="Ganglioneuroma"),"y"]-4
  lab[which(lab$Diagnosis=="Ganglioneuroma"),"x"]<-lab[which(lab$Diagnosis=="Ganglioneuroma"),"x"]-4
  lab[which(lab$Diagnosis=="SFT"),"x"]<-lab[which(lab$Diagnosis=="SFT"),"x"]+4.5
  lab[which(lab$Diagnosis=="SFT"),"y"]<-lab[which(lab$Diagnosis=="SFT"),"y"]-2
  lab[which(lab$Diagnosis=="PH/PG"),"x"]<-lab[which(lab$Diagnosis=="PH/PG"),"x"]+12
  lab[which(lab$Diagnosis=="Pituitary"),"y"]<-lab[which(lab$Diagnosis=="Pituitary"),"y"]+3.5
  lab[which(lab$Diagnosis=="PitNET"),"y"]<-lab[which(lab$Diagnosis=="PitNET"),"y"]-3.5
  lab[which(lab$Diagnosis=="MB-GP3"),"x"]<-lab[which(lab$Diagnosis=="MB-GP3"),"x"]+13
  lab[which(lab$Diagnosis=="MB-GP4"),"x"]<-lab[which(lab$Diagnosis=="MB-GP4"),"x"]-14
  lab[which(lab$Diagnosis=="MB-WNT"),"y"]<-lab[which(lab$Diagnosis=="MB-WNT"),"y"]-5
  lab[which(lab$Diagnosis=="MB-SHH"),"x"]<-lab[which(lab$Diagnosis=="MB-SHH"),"x"]+14
  lab[which(lab$Diagnosis=="MB-SHH"),"y"]<-lab[which(lab$Diagnosis=="MB-SHH"),"y"]-4
  lab[which(lab$Diagnosis=="Peri. Nerve"),"y"]<-lab[which(lab$Diagnosis=="Peri. Nerve"),"y"]+4
  lab[which(lab$Diagnosis=="HGNET-MN1"),"y"]<-lab[which(lab$Diagnosis=="HGNET-MN1"),"y"]+4
  lab[which(lab$Diagnosis=="Retina, Fetal"),"y"]<-lab[which(lab$Diagnosis=="Retina, Fetal"),"y"]-2.5
  lab[which(lab$Diagnosis=="Retina, Fetal"),"x"]<-lab[which(lab$Diagnosis=="Retina, Fetal"),"x"]-10
  lab[which(lab$Diagnosis=="Retina"),"x"]<-lab[which(lab$Diagnosis=="Retina"),"x"]-7
  lab[which(lab$Diagnosis=="Retinoblastoma"),"y"]<-lab[which(lab$Diagnosis=="Retinoblastoma"),"y"]-8
  lab[which(lab$Diagnosis=="Nerve Ganglia"),"y"]<-lab[which(lab$Diagnosis=="Nerve Ganglia"),"y"]+5
  lab[which(lab$Diagnosis=="Cerebellum"),"y"]<-lab[which(lab$Diagnosis=="Cerebellum"),"y"]-4
  lab[which(lab$Diagnosis=="Cerebellum, Fetal"),"y"]<-lab[which(lab$Diagnosis=="Cerebellum, Fetal"),"y"]+3
  lab[which(lab$Diagnosis=="Cerebellum, Fetal"),"x"]<-lab[which(lab$Diagnosis=="Cerebellum, Fetal"),"x"]-5
  lab[which(lab$Diagnosis=="NB-FOXR2"),"y"]<-lab[which(lab$Diagnosis=="NB-FOXR2"),"y"]-4
  lab[which(lab$Diagnosis=="ETMR"),"x"]<-lab[which(lab$Diagnosis=="ETMR"),"x"]+7
  lab[which(lab$Diagnosis=="HGNET-BCOR"),"y"]<-lab[which(lab$Diagnosis=="HGNET-BCOR"),"y"]-3.5
  lab[which(lab$Diagnosis=="HGNET-BCOR"),"x"]<-lab[which(lab$Diagnosis=="HGNET-BCOR"),"x"]-5
  lab[which(lab$Diagnosis=="IDH Mut-1p19q del"),"y"]<-lab[which(lab$Diagnosis=="IDH Mut-1p19q del"),"y"]-10
  lab[which(lab$Diagnosis=="IDH Mut-1p19q del"),"x"]<-lab[which(lab$Diagnosis=="IDH Mut-1p19q del"),"x"]-17.5
  lab[which(lab$Diagnosis=="AT/RT-SHH"),"x"]<-lab[which(lab$Diagnosis=="AT/RT-SHH"),"x"]+11
  lab[which(lab$Diagnosis=="PCNSL"),"x"]<-lab[which(lab$Diagnosis=="PCNSL"),"x"]+8.5
  lab[which(lab$Diagnosis=="PCNSL"),"y"]<-lab[which(lab$Diagnosis=="PCNSL"),"y"]+2
  lab[which(lab$Diagnosis=="Midline"),"x"]<-lab[which(lab$Diagnosis=="Midline"),"x"]+11.5
  lab[which(lab$Diagnosis=="Midline"),"y"]<-lab[which(lab$Diagnosis=="Midline"),"y"]-8
  lab[which(lab$Diagnosis=="IDH Mutant"),"x"]<-lab[which(lab$Diagnosis=="IDH Mutant"),"x"]+12.75
  lab[which(lab$Diagnosis=="IDH Mutant"),"y"]<-lab[which(lab$Diagnosis=="IDH Mutant"),"y"]-16
  lab[which(lab$Diagnosis=="G34 Mutant"),"x"]<-lab[which(lab$Diagnosis=="G34 Mutant"),"x"]+10
  lab[which(lab$Diagnosis=="G34 Mutant"),"y"]<-lab[which(lab$Diagnosis=="G34 Mutant"),"y"]-3
  lab[which(lab$Diagnosis=="Lung Met"),"x"]<-lab[which(lab$Diagnosis=="Lung Met"),"x"]-6
  lab[which(lab$Diagnosis=="Lung Met"),"y"]<-lab[which(lab$Diagnosis=="Lung Met"),"y"]-3.5
  lab[which(lab$Diagnosis=="Breast Met"),"x"]<-lab[which(lab$Diagnosis=="Breast Met"),"x"]+9
  lab[which(lab$Diagnosis=="Ada. Cranio."),"x"]<-lab[which(lab$Diagnosis=="Ada. Cranio."),"x"]-8
  lab[which(lab$Diagnosis=="Ada. Cranio."),"y"]<-lab[which(lab$Diagnosis=="Ada. Cranio."),"y"]+2.5
  lab[which(lab$Diagnosis=="GCT"),"x"]<-lab[which(lab$Diagnosis=="GCT"),"x"]-4
  lab[which(lab$Diagnosis=="Neurofibroma"),"x"]<-lab[which(lab$Diagnosis=="Neurofibroma"),"x"]+8
  lab[which(lab$Diagnosis=="Neurofibroma"),"y"]<-lab[which(lab$Diagnosis=="Neurofibroma"),"y"]+4
  lab[which(lab$Diagnosis=="MPNST"),"y"]<-lab[which(lab$Diagnosis=="MPNST"),"y"]+2.25
  lab[which(lab$Diagnosis=="MPNST"),"x"]<-lab[which(lab$Diagnosis=="MPNST"),"x"]+6.5
  lab[which(lab$Diagnosis=="EFT-CIC"),"y"]<-lab[which(lab$Diagnosis=="EFT-CIC"),"y"]+1
  lab[which(lab$Diagnosis=="EFT-CIC"),"x"]<-lab[which(lab$Diagnosis=="EFT-CIC"),"x"]+8
  lab[which(lab$Diagnosis=="AT/RT-TYR"),"x"]<-lab[which(lab$Diagnosis=="AT/RT-TYR"),"x"]+10
  lab[which(lab$Diagnosis=="AT/RT-TYR"),"y"]<-lab[which(lab$Diagnosis=="AT/RT-TYR"),"y"]+2.5
  lab[which(lab$Diagnosis=="AT/RT-MYC"),"x"]<-lab[which(lab$Diagnosis=="AT/RT-MYC"),"x"]+4
  lab[which(lab$Diagnosis=="AT/RT-MYC"),"y"]<-lab[which(lab$Diagnosis=="AT/RT-MYC"),"y"]-6
  lab[which(lab$Diagnosis=="PXA"),"y"]<-lab[which(lab$Diagnosis=="PXA"),"y"]+6
  lab[which(lab$Diagnosis=="PXA"),"x"]<-lab[which(lab$Diagnosis=="PXA"),"x"]+3
  lab[which(lab$Diagnosis=="Angio."),"x"]<-lab[which(lab$Diagnosis=="Angio."),"x"]-5
  lab[which(lab$Diagnosis=="Angio."),"y"]<-lab[which(lab$Diagnosis=="Angio."),"y"]+12
  lab[which(lab$Diagnosis=="SEGA"),"x"]<-lab[which(lab$Diagnosis=="SEGA"),"x"]-12
  lab[which(lab$Diagnosis=="SEGA"),"y"]<-lab[which(lab$Diagnosis=="SEGA"),"y"]-4
  lab[which(lab$Diagnosis=="DIG"),"x"]<-lab[which(lab$Diagnosis=="DIG"),"x"]+4
  lab[which(lab$Diagnosis=="DIG"),"y"]<-lab[which(lab$Diagnosis=="DIG"),"y"]+8
  lab[which(lab$Diagnosis=="GG"),"x"]<-lab[which(lab$Diagnosis=="GG"),"x"]-9
  lab[which(lab$Diagnosis=="GG"),"y"]<-lab[which(lab$Diagnosis=="GG"),"y"]-2
  lab[which(lab$Diagnosis=="PA"),"x"]<-lab[which(lab$Diagnosis=="PA"),"x"]-8
  lab[which(lab$Diagnosis=="PA"),"y"]<-lab[which(lab$Diagnosis=="PA"),"y"]+7
  lab[which(lab$Diagnosis=="PFA"),"x"]<-lab[which(lab$Diagnosis=="PFA"),"x"]-8
  lab[which(lab$Diagnosis=="PFA"),"y"]<-lab[which(lab$Diagnosis=="PFA"),"y"]-8
  lab[which(lab$Diagnosis=="Meningioma"),"y"]<-lab[which(lab$Diagnosis=="Meningioma"),"y"]+9
  lab[which(lab$Diagnosis=="Uveal Melanoma"),"y"]<-lab[which(lab$Diagnosis=="Uveal Melanoma"),"y"]+7
}

{
  p<-ggplot()+geom_point(data=d, aes(tsne1,tsne2,color=factor(Diagnosis)), size=0.75)+theme_pubr()+
    geom_text(data = lab, aes(x=x,y=y,label=Diagnosis), size=3)+
    annotate("text", x=-94, y=82, label= "t-SNE Dimensionality Reduction\nAfter Classification\n(7334 samples)", hjust = 0, fontface = "bold") + 
    annotate("text", x=-94, y=-75, size=3, label= "Angio. = Angiocentric\nGCT = Germ Cell Tumor\nGG = Ganglioglioma\nSEGA = Subependymal Giant\n              Cell Astrocytoma", hjust = 0) + 
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

#edited in illustrator 