library(ggplot2)
library(ggnewscale)
library(cowplot)
library(dplyr)
library(ggpubr)
library(plotly)

d<-read.csv("Atlas_data.csv", header=T, row.names = 1, stringsAsFactors = T)

d<-d[d$Final!="-",]

d$Final<-factor(d$Final, levels = c("Papilloma", "Angiocentric", "Diffuse Glioma", "G34 Mutant", "IDH Mutant", "Midline",
                                    "HGNET-MN1", "PA", "PXA", "AT/RT-MYC", "AT/RT-SHH", "AT/RT-TYR", "ETMR", "HGNET-BCOR","MB-GP3",
                                    "MB-GP4","MB-SHH","MB-WNT", "Neuroblastoma", "NB-FOXR2", "Retinoblastoma", "MPE", "PFA", "PFB", "RELA", "SE", "YAP", "Germ Cell Tumor",
                                    "DIG", "Ganglioglioma", "PCNSL", "Uveal Melanoma", "EFT-CIC", "Meningioma", "SFT", "Ganglioneuroma",
                                    "MPNST", "Neurofibroma", "PH/PG", "Ada. Cranio.", "PitNET", "Breast", "Lung",
                                    "Cerebellum", "Cerebellum, Fetal", "Choroid Plexus", "CNS", "Nerve Ganglia", "Peri. Nerve",
                                    "Pituitary", "Retina", "Retina, Fetal"))
#https://medialab.github.io/iwanthue/ #44 colors
colors<-c("#6a0033",
          "#e4016e",
          "#ff527b",
          "#360e15",
          "#ff757c",
          "#dc002d",
          "#ff9170",
          "#7a2000",
          "#ff7e40",
          "#ac5300",
          "#da6c00",
          "#401f00",
          "#ffb173",
          "#d79600",
          "#6a4800",
          "#d6c6b2",
          "#ddaf00",
          "#d6ca6f",
          "#6e6900",
          "#1d1c10",
          "#8db600",
          "#1b7a00",
          "#8ddb76",
          "#00b861",
          "#9ad5b1",
          "#005d3f",
          "#01bbb7",
          "#00444d",
          "#01afc7",
          "#54d8f9",
          "#0189dd",
          "#8ab5ff",
          "#5292ff",
          "#004690",
          "#00317f",
          "#4263ff",
          "#240a4e",
          "#271332",
          "#fa63ff",
          "#760078",
          "#ff77f6",
          "#8a005f",
          "#ffa5ca",
          "black","grey10","grey20","gray80","grey40","grey50","grey60","grey70","grey30") #can add one more "#ff4494"

lab<- d %>% group_by(Final) %>% summarize(x = mean(V1), y = mean(V2))

lab[1,2]<-lab[1,2]-7 #papilloma x
lab[1,3]<-lab[1,3]-3 #papilloma y
lab[2,2]<-lab[2,2]   #angio x
lab[2,3]<-lab[2,3]+3 #angio y
lab[3,2]<-lab[3,2]-27   #diffuse glioma x
lab[3,3]<-lab[3,3]-22   #diffuse glioma y
lab[4,2]<-lab[4,2]+3   #G34 x
lab[4,3]<-lab[4,3]+5      #G34 y
lab[5,2]<-lab[5,2]-8   #IDH x
lab[5,3]<-lab[5,3]+16   #IDH y
lab[6,2]<-lab[6,2]-15   #Midline x
lab[6,3]<-lab[6,3]+15   #Midline y
lab[7,2]<-lab[7,2]   #MN1 x
lab[7,3]<-lab[7,3]-4   #MN1 y
lab[8,2]<-lab[8,2]+13   #PA x
lab[8,3]<-lab[8,3]+6   #PA y
lab[9,2]<-lab[9,2]+4   #PXA x
lab[9,3]<-lab[9,3]+1   #PXA y
lab[10,2]<-lab[10,2]+12   #AT MYC x
lab[10,3]<-lab[10,3]-2   #AT MYC y
lab[11,2]<-lab[11,2]+12   #AT SHH x
lab[11,3]<-lab[11,3]+3   #AT SHH y
lab[12,2]<-lab[12,2]+13   #AT TYR x
lab[12,3]<-lab[12,3]+2   #AT TYR y
lab[13,2]<-lab[13,2]+6   #ETMR x
lab[13,3]<-lab[13,3]+5   #ETMR y
lab[14,2]<-lab[14,2]+13   #BCOR x
lab[14,3]<-lab[14,3]-1   #BCOR y
lab[15,2]<-lab[15,2]+12   #MB3 x
lab[15,3]<-lab[15,3]+7   #MB3 y
lab[16,2]<-lab[16,2]-15   #MB4 x
lab[16,3]<-lab[16,3]+10   #MB4 y
lab[17,2]<-lab[17,2]+16   #MBSHH x
lab[17,3]<-lab[17,3]   #MBSHH y
lab[18,2]<-lab[18,2]+12   #MBWNT x
lab[18,3]<-lab[18,3]   #MBWNT y
lab[19,2]<-lab[19,2]  #Neuroblastoma x
lab[19,3]<-lab[19,3]-10   #Neuroblastoma y
lab[20,2]<-lab[20,2]   #NB-FOXR2 x
lab[20,3]<-lab[20,3]+4   #NB-FOXR2 y
lab[21,2]<-lab[21,2]-20   #Retinoblastoma x
lab[21,3]<-lab[21,3]   #Retinoblastoma y
lab[22,2]<-lab[22,2]+8   #MPE x
lab[22,3]<-lab[22,3]   #MPE y
lab[23,2]<-lab[23,2]-6   #PFA x
lab[23,3]<-lab[23,3]-10   #PFA y
lab[24,2]<-lab[24,2]+6   #PFB x
lab[24,3]<-lab[24,3]-4   #PFB y
lab[25,2]<-lab[25,2]+11   #RELA x
lab[25,3]<-lab[25,3]   #RELA y
lab[26,2]<-lab[26,2]+3   #SE x
lab[26,3]<-lab[26,3]+6   #SE y
lab[27,2]<-lab[27,2]+6   #YAP x
lab[27,3]<-lab[27,3]-2   #YAP y
lab$Final<-as.character(lab$Final)
lab[2,1]<-"Angio."
lab[28,1]<-"GCT"
lab[32,1]<-"UM"
lab[36,1]<-"GN"
lab[4,1]<-"G34"
lab[42,1]<-"BrMet"
lab[43,1]<-"LungMet"
lab[44,1]<-"Cerbllm"
lab$Final<-factor(lab$Final, levels=unique(lab$Final))
lab[28,2]<-lab[28,2]-4   #GCT x
lab[28,3]<-lab[28,3]   #GCT y
lab[29,2]<-lab[29,2]+8.5   #DIG x
lab[29,3]<-lab[29,3]-3   #DIG y
lab[30,2]<-lab[30,2]   #GG x
lab[30,3]<-lab[30,3]+5   #GG y
lab[31,2]<-lab[31,2]+11   #PCNSL x
lab[31,3]<-lab[31,3]   #PCNSL y
lab[32,2]<-lab[32,2]+10   #UM x
lab[32,3]<-lab[32,3]   #UM y
lab[33,2]<-lab[33,2]+6   #CIC x
lab[33,3]<-lab[33,3]-2   #CIC y
lab[34,2]<-lab[34,2]+10   #Mening x
lab[34,3]<-lab[34,3]-10   #Mening y
lab[35,2]<-lab[35,2]+5   #SFT x
lab[35,3]<-lab[35,3]   #SFT y
lab[36,2]<-lab[36,2]-2   #ganglioneuroma x
lab[36,3]<-lab[36,3]+5   #ganglioneuroma y
lab[37,2]<-lab[37,2]+7   #MPNST x
lab[37,3]<-lab[37,3]-4   #MPSNT y
lab[38,2]<-lab[38,2]+14   #NF x
lab[38,3]<-lab[38,3]   #NF y
lab[39,2]<-lab[39,2]+13   #PHPG x
lab[39,3]<-lab[39,3]   #PHPG y
lab[40,2]<-lab[40,2]   #Cranio x
lab[40,3]<-lab[40,3]-3   #Crnio y
lab[41,2]<-lab[41,2]+8   #PitNET x
lab[41,3]<-lab[41,3]+3   #PitNET y
lab[42,2]<-lab[42,2]-7   #Breast x
lab[42,3]<-lab[42,3]-5   #Breast y
lab[43,2]<-lab[43,2]+2   #Lung x main
lab[43,3]<-lab[43,3]-5   #Lung y main
lab[44,2]<-lab[44,2]-18   #cerebellum x
lab[44,3]<-lab[44,3]-1   #cerebellum y
lab[45,2]<-lab[45,2]-14   #cerebellum,fetal x
lab[45,3]<-lab[45,3]   #cerebellum,fetal y
lab[46,2]<-lab[46,2]-10   #Ch.plex x
lab[46,3]<-lab[46,3]+3   #Ch.plex y
lab[47,2]<-lab[47,2]   #CNS x
lab[47,3]<-lab[47,3]+40   #CNS y
lab[48,2]<-lab[48,2]   #Ganglia x
lab[48,3]<-lab[48,3]-5   #Ganglia y
lab[49,2]<-lab[49,2]+8   #PN x
lab[49,3]<-lab[49,3]-4   #PN y
lab[50,2]<-lab[50,2]+73   #pit x
lab[50,3]<-lab[50,3]+32   #pit y
lab[51,2]<-lab[51,2]   #Retina x
lab[51,3]<-lab[51,3]-4   #Retina y
lab[52,2]<-lab[52,2]+9   #Fetal Retina x
lab[52,3]<-lab[52,3]-3   #Fetal Retina y
lab<-rbind(lab,c("Pituitary",-90,-38))
lab<-rbind(lab,c("LungMet",56,-27))
lab<-rbind(lab,c("Cerbllm",-43,1))
lab$Final<-factor(lab$Final, levels=unique(lab$Final))
lab$x<-as.numeric(lab$x)
lab$y<-as.numeric(lab$y)

p<-ggplot()+geom_point(data=d, aes(V1,V2,color=factor(Final)), size=0.75)+theme_pubr()+
  geom_text(data = lab, aes(x=x,y=y,label=Final), size=3)+
  annotate("text", x=-94, y=87, label= "t-SNE Dimensionality Reduction\n(7296 samples)", hjust = 0, fontface = "bold") + 
  annotate("text", x=-94, y=-75, size=3, label= "Angio. = Angiocentric\nBrMet = Breast Cancer Metastasis\nCerbllm = Cerebellum\nGCT = Germ Cell Tumor\nGN = Ganglioneuroma\nUM = Uveal Melanoma", hjust = 0) + 
  scale_color_manual(values=colors)+
  theme(legend.position = "none")+rremove("axis")+rremove("axis.text")+rremove("ticks")+rremove("xylab")+
  theme(panel.border = element_rect(colour = "black", fill=NA))
#p

p1<-ggplot()+
  geom_point(data=d[d$Class2=="Choroid Plexus",],aes(V1,V2,color=Final),size=1)+
  scale_color_manual(values=colors[1], name="Choroid Plexus")+
  guides(color = guide_legend(override.aes=list(shape = 15, size=5), order = 1))+
  new_scale_color() +
  
  geom_point(data=d[d$Class2=="Diffuse Glioma",],aes(V1,V2,color=Final),size=1)+
  scale_color_manual(values=colors[2:6], name="Diffuse Glioma")+
  guides(color = guide_legend(override.aes=list(shape = 15, size=5, ncol=2), order = 2))+
  new_scale_color() +
  
  geom_point(data=d[d$Class2=="Discrete Glioma",],aes(V1,V2,color=Final),size=1)+
  scale_color_manual(values=colors[7:9], name="Discrete Glioma")+
  guides(color = guide_legend(override.aes=list(shape = 15, size=5), order = 3))+
  new_scale_color() +
  
  geom_point(data=d[d$Class2=="Embryonal",],aes(V1,V2,color=Final),size=1)+
  scale_color_manual(values=colors[10:21], name="Embryonal")+
  guides(color = guide_legend(override.aes=list(shape = 15, size=5), order = 4))+
  
  theme_classic()+theme(legend.title=element_text(face="bold"), legend.justification = "top")

p1l <- get_legend(p1)

p2<-ggplot()+
  geom_point(data=d[d$Class2=="Ependymal",],aes(V1,V2,color=Final),size=1)+
  scale_color_manual(values=colors[22:27], name="Ependymal")+
  guides(color = guide_legend(override.aes=list(shape = 15, size=5), order = 1))+
  new_scale_color() +
  
  geom_point(data=d[d$Class2=="Germ Cell",],aes(V1,V2,color=Final),size=1)+
  scale_color_manual(values=colors[28], name="Germ Cell")+
  guides(color = guide_legend(override.aes=list(shape = 15, size=5), order = 2))+
  new_scale_color() +
  
  geom_point(data=d[d$Class2=="Glioneuronal",],aes(V1,V2,color=Final),size=1)+
  scale_color_manual(values=colors[29:30], name="Glioneuronal")+
  guides(color = guide_legend(override.aes=list(shape = 15, size=5), order = 3))+
  new_scale_color() +
  
  geom_point(data=d[d$Class2=="Hematolymphoid",],aes(V1,V2,color=Final),size=1)+
  scale_color_manual(values=colors[31], name="Hematolymphoid")+
  guides(color = guide_legend(override.aes=list(shape = 15, size=5), order = 4))+
  new_scale_color() +
  
  geom_point(data=d[d$Class2=="Melanocytic",],aes(V1,V2,color=Final),size=1)+
  scale_color_manual(values=colors[32], name="Melanocytic")+
  guides(color = guide_legend(override.aes=list(shape = 15, size=5), order = 5))+
  new_scale_color() +
  
  geom_point(data=d[d$Class2=="Mesenchymal",],aes(V1,V2,color=Final),size=1)+
  scale_color_manual(values=colors[33:35], name="Mesenchymal")+
  guides(color = guide_legend(override.aes=list(shape = 15, size=5), order = 6))+
  
  theme_classic()+theme(legend.title=element_text(face="bold"),legend.justification = "top")

p2l <- get_legend(p2)

p3<-ggplot()+
  geom_point(data=d[d$Class2=="Nerve",],aes(V1,V2,color=Final),size=1)+
  scale_color_manual(values=colors[36:39], name="Nerve")+
  guides(color = guide_legend(override.aes=list(shape = 15, size=5), order = 1))+
  new_scale_color() +
  
  geom_point(data=d[d$Class2=="Sellar",],aes(V1,V2,color=Final),size=1)+
  scale_color_manual(values=colors[40:41], name="Sellar")+
  guides(color = guide_legend(override.aes=list(shape = 15, size=5), order = 2))+
  new_scale_color() +
  
  geom_point(data=d[d$Class2=="Metastasis",],aes(V1,V2,color=Final),size=1)+
  scale_color_manual(values=colors[42:46], name="Metastasis")+
  guides(color = guide_legend(override.aes=list(shape = 15, size=5), order = 3))+
  new_scale_color() +
  
  geom_point(data=d[d$Class2=="Non-tumor",],aes(V1,V2,color=Final),size=1)+
  scale_color_manual(values=c("black","grey10","grey20","gray80","grey40","grey50","grey60","grey70","grey30"), name="Non-tumor")+
  guides(color = guide_legend(override.aes=list(shape = 15, size=5), order = 4))+
  new_scale_color() +
  
  theme_classic()+theme(legend.title=element_text(face="bold"),legend.justification = "top")

p3l <- get_legend(p3)

plot_grid(p,p1l,p2l,p3l, ncol=4, rel_widths = c(6,1,1,1), align = 'hv')
#1200 x 670 svg
