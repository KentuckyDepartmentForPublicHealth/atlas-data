library(ggpubr)
library(ggplot2)
library(cowplot)
library(forcats)

d<-read.csv("~/Library/CloudStorage/OneDrive-uoflhealth/Atlas/Atlas_data.csv", header=T, stringsAsFactors = T)
d$Compartment<-factor(d$Compartment, levels = c("Brain, Supratentorial","Brain, Infratentorial","Peripheral","Eye","Spinal cord","Brain","NA"), 
                      labels = c("Supratentorial","Infratentorial","Peripheral","Eye","Spinal cord","Brain, NOS","NA"))
d$Age..grouped.<-factor(d$Age..grouped., levels = c("Fetal","0-5yrs","5-10yrs","10-20yrs","20-40yrs","40-60yrs","60-80yrs","80+yrs","NA"),
                        labels = c("Fetal","0-5 yrs","5-10 yrs","10-20 yrs","20-40 yrs","40-60 yrs","60-80 yrs","80+ yrs","NA"))

d$Country..listed.on.deposited.data.<-fct_infreq(d$Country..listed.on.deposited.data.)

#ggplot(d)+geom_bar(aes(x=reorder(Country..listed.on.deposited.data.,Country..listed.on.deposited.data.,function(x)-length(x))))+
#          geom_bar(aes(x=reorder(Compartment,Compartment,function(x)-length(x))))+
#          geom_bar(aes(x=Age..grouped., fill=Sex))+scale_fill_manual(values=c("hotpink", "blue"))+
#          theme_pubr()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab("Number of Samples")


library(dplyr)
comp<-d %>% group_by(Compartment,Sex) %>% tally()
colnames(comp)[1]<-c("var")
comp$type<-"Location"

cntry<-d %>% group_by(Country..listed.on.deposited.data.) %>% tally()
colnames(cntry)[1]<-c("var")
cntry$type<-"Geographic Representation"
cntry$Sex<-NA

age<-d %>% group_by(Age..grouped.,Sex) %>% tally()
colnames(age)[1]<-c("var")
age$type<-"Age"

library(tidyr)
d1<-d[,c(40:47)]
d1<-gather(d1, gene, mut, X1p.19q.codel:MCYN.amplification, factor_key=TRUE)
d1$mut<-factor(d1$mut)

d1<-subset(d1, mut!="NA")
d1<-d1 %>% group_by(gene,mut) %>% tally()
d1<-d1[-c(11,18),]
d1$mut<-gsub("G34","Yes",d1$mut)
d1$mut<-gsub("K27","Yes",d1$mut)
d1$mut<-gsub("MUT","Yes",d1$mut)
d1$mut<-gsub("WT","No",d1$mut)
d1$mut<-gsub("Methylated","Yes",d1$mut)
d1$mut<-gsub("Unmethylated","No",d1$mut)
d1$mut<-gsub("KIAA1549-BRAF","Yes",d1$mut)
d1$type<-"Genetic Profile"

d1$gene<-factor(d1$gene, levels=c("X1p.19q.codel","IDH1.2.mutation","H3.mutation", "TERT..promoter..mutation", "EGFR.amplification",
                                  "MGMT.promoter.methylation", "BRAF.mutation", "MCYN.amplification"),
                labels=c("1p/19q codel.","IDH1/2 mut.","H3 mut.", "pTERT mut.", "EGFR amplif.",
                         "pMGMT meth.", "BRAF mut.", "MCYN amplif."))

colnames(d1)[1:2]<-c("var","Sex")

newd<-rbind(cntry,comp,age,d1)
newd$type<-factor(newd$type, levels=c("Geographic Representation","Age","Location","Genetic Profile"))
newd$Sex<-gsub("F","Female",newd$Sex)
newd$Sex<-gsub("M","Male",newd$Sex)

A<-ggplot(newd, aes(x = var,  y = n, fill = Sex))+scale_fill_manual(values=c("hotpink","blue","darkgreen","red","grey"))+
  geom_hline(yintercept=250, linetype="dashed", color = "grey")+
  geom_hline(yintercept=500, linetype="dashed", color = "grey")+
  geom_col(position = "stack")+scale_y_continuous(breaks=c(0, 250, 500, 1000, 2000, 3000))+
  facet_grid(~type, scales = "free", space = "free")+theme_pubr()+xlab("")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.title=element_blank(),panel.spacing = unit(1.5, "lines"),
        legend.position = c(0.75, 0.9),legend.direction = "horizontal", strip.text = element_text(size=14),
        axis.title=element_text(size=14))+
  ylab("Number of Samples")

library(survival)
library(survminer)

d<-read.csv("Atlas_data.csv", header=T)
d<-d[grep("^MB",d$Final.diagnosis_bc),]
fit<- survfit(Surv(as.numeric(Overall.survival..months.), Vital.status..1.dead..0.alive.) ~ Final.diagnosis_bc, data = d)
mb<-ggsurvplot(fit, risk.table = F, legend='none', xlab="Time (Months)")$plot+
  annotate("text", x=78, y=.95, label= "MB-WNT (n=18)",hjust = 0)+
  annotate("text", x=78, y=.79, label= "MB-GP4 (n=48)",hjust = 0)+
  annotate("text", x=78, y=.70, label= "MB-GP3 (n=33)",hjust = 0)+
  annotate("text", x=78, y=.55, label= "MB-SHH (n=24)",hjust = 0)

d<-read.csv("Atlas_data.csv", header=T)
d<-d[grep("^AT",d$Final.diagnosis_bc),]
fit<- survfit(Surv(as.numeric(Overall.survival..months.), Vital.status..1.dead..0.alive.) ~ Final.diagnosis_bc, data = d)
at<-ggsurvplot(fit, risk.table = F, legend='none', xlab="Time (Months)",xlim=c(0,50),break.x.by=10)$plot+
  annotate("text", x=25, y=.82, label= "AT/RT-MYC (n=4)",hjust = 0)+
  annotate("text", x=25, y=.34, label= "AT/RT-SHH (n=10)",hjust = 0)+
  annotate("text", x=25, y=.08, label= "AT/RT-TYR (n=9)",hjust = 0)

d<-read.csv("Atlas_data.csv", header=T)
d<-d[d$Final.diagnosis_bc=="ETMR"|d$Final.diagnosis_bc=="HGNET-BCOR"|d$Final.diagnosis_bc=="NB-FOXR2",]
fit<- survfit(Surv(as.numeric(Overall.survival..months.), Vital.status..1.dead..0.alive.) ~ Final.diagnosis_bc, data = d)
em<-ggsurvplot(fit, risk.table = F, legend='none', xlab="Time (Months)")$plot+
  annotate("text", x=35, y=.93, label= "NB-FOXR2 (n=6)",hjust = 0)+
  annotate("text", x=35, y=.56, label= "HGNET-BCOR (n=9)",hjust = 0)+
  annotate("text", x=35, y=.17, label= "ETMR (n=10)",hjust = 0)

d<-read.csv("Atlas_data.csv", header=T)
d<-d[d$Final.diagnosis_bc=="Neuroblastoma",]
fit<- survfit(Surv(as.numeric(Overall.survival..months.), Vital.status..1.dead..0.alive.) ~ MCYN.amplification, data = d)
nb<-ggsurvplot(fit, risk.table = F, legend='none', xlab="Time (Months)")$plot+
  annotate("text", x=50, y=.85, label= "NB: MCYN Not Amplified (n=71)",hjust = 0)+
  annotate("text", x=55, y=.20, label= "NB: MCYN Amplification (n=16)",hjust = 0)

d<-read.csv("Atlas_data.csv", header=T)
d<-d[d$Final.diagnosis_bc=="YAP"|d$Final.diagnosis_bc=="PFB"|d$Final.diagnosis_bc=="PFA"|d$Final.diagnosis_bc=="RELA",]
fit<- survfit(Surv(as.numeric(Overall.survival..months.), Vital.status..1.dead..0.alive.) ~ Final.diagnosis_bc, data = d)
ep<-ggsurvplot(fit, risk.table = T, legend='none', xlab="Time (Months)")$plot+
  annotate("text", x=90, y=1, label= "PFB (n=18)/YAP (n=7)",hjust = 0)+
  annotate("text", x=105, y=.51, label= "PFA (n=148)",hjust = 0)+
  annotate("text", x=93, y=.22, label= "RELA (n=18)",hjust = 0)

d<-read.csv("Atlas_data.csv", header=T)
d<-d[d$Final.diagnosis_bc=="Meningioma",]
d$Grade<-gsub("4","3",d$Grade)
fit<- survfit(Surv(as.numeric(Overall.survival..months.), Vital.status..1.dead..0.alive.) ~ Grade, data = d)
men<-ggsurvplot(fit, risk.table = F, legend='none', xlab="Time (Months)")$plot+
  annotate("text", x=55, y=1, label= "MEN Grade 1 (n=43)",hjust = 0)+
  annotate("text", x=37, y=.75, label= "MEN Grade 2\n(n=18)",hjust = 0)+
  annotate("text", x=17, y=.43, label= "MEN Grade 3\n(n=8)",hjust = 0)

d<-read.csv("Atlas_data.csv", header=T)
d<-d[d$Final.diagnosis_bc=="Diffuse Glioma-WT"|d$Final.diagnosis_bc=="IDH Mutant"|d$Final.diagnosis_bc=="IDH Mut-1p19q del"|d$Final.diagnosis_bc=="G34 Mutant"|d$Final.diagnosis_bc=="Midline",]
d<-d[d$Dataset.ID!="GSE53733",]
fit<- survfit(Surv(as.numeric(Overall.survival..months.), Vital.status..1.dead..0.alive.) ~ Final.diagnosis_bc, data = d)
dg<-ggsurvplot(fit, risk.table = F, legend='none', xlab="Time (Months)")$plot+
  annotate("text", x=130, y=.8, label= "IDH mutant (n=237)",hjust = 0, vjust = 0, color ="#00B0F6" )+
  annotate("text", x=130, y=.9, label= "IDH Mut-1p19q del\n             (n=284)",hjust = 0,vjust = 0,color="#00BF7D" )+
  annotate("text", x=130, y=.65, label= "G34 mutant (n=15)",hjust = 0, vjust = 0,color="#A3A500")+
  annotate("text", x=130, y=.3, label= "Diffuse glioma-WT\n             (n=795)",hjust = 0,vjust = 0, color="#F8766D")+
  annotate("text", x=130, y=.5, label= "Midline (n=57)",hjust = 0, vjust = 0,color="#E76BF3")

d<-read.csv("Atlas_data.csv", header=T)
d<-d[d$Final.diagnosis_bc=="Diffuse Glioma-WT"|d$Final.diagnosis_bc=="IDH Mutant"|d$Final.diagnosis_bc=="IDH Mut-1p19q del"|d$Final.diagnosis_bc=="G34 Mutant"|d$Final.diagnosis_bc=="Midline",]
d<-d[d$Dataset.ID!="GSE53733",]
d<-d[d$Grade!="Grade 1",]
fit<- survfit(Surv(as.numeric(Overall.survival..months.), Vital.status..1.dead..0.alive.) ~ Grade, data = d)
dgg<-ggsurvplot(fit, risk.table = F, legend='none', xlab="Time (Months)")$plot+
  annotate("text", x=55, y=.75, label= "Grade 2 (n=239)",hjust = 0)+
  annotate("text", x=145, y=.25, label= "Grade 3 (n=411)",hjust = 0,color="#00BA38")+
  annotate("text", x=50, y=.2, label= "Grade 4\n(n=672)",hjust = 0)

d<-read.csv("Atlas_data.csv", header=T)
d<-d[d$Final.diagnosis_bc=="Diffuse Glioma-WT"|d$Final.diagnosis_bc=="IDH Mutant"|d$Final.diagnosis_bc=="IDH Mut-1p19q del"|d$Final.diagnosis_bc=="G34 Mutant"|d$Final.diagnosis_bc=="Midline",]
d<-d[d$Dataset.ID!="GSE53733",]
d<-d[d$Grade!="Grade 1",]
fit<- survfit(Surv(as.numeric(Overall.survival..months.), Vital.status..1.dead..0.alive.) ~ MGMT.promoter.methylation, data = d)
mgmt<-ggsurvplot(fit, risk.table = F, legend='none', xlab="Time (Months)")$plot+
  annotate("text", x=50, y=.75, label= "MGMT Promoter Methylated\n(n=175)",hjust = 0)+
  annotate("text", x=47, y=.18, label= "MGMT Promoter Unmethylated\n                                   (n=78)",hjust = 0)

d<-read.csv("Atlas_data.csv", header=T)
d<-d[d$Final.diagnosis_bc=="PA"|d$Final.diagnosis_bc=="Germ Cell Tumor",]
fit<- survfit(Surv(as.numeric(Overall.survival..months.), Vital.status..1.dead..0.alive.) ~ Final.diagnosis_bc, data = d)
pa<-ggsurvplot(fit, risk.table = F, legend='none', xlab="Time (Months)")$plot+
  annotate("text", x=45, y=.75, label= "PA (n=28)",hjust = 0)+
  annotate("text", x=45, y=1, label= "GCT (n=11)",hjust = 0)
#350x300 each plot size above

library(cowplot)
B<-plot_grid(dg,dgg,mgmt,mb,at,nb,em,ep,men,pa, nrow = 2, labels = c("E","F","G","H","I","J","K","L","M","N"), align = 'hv',label_size = 16)
#1750x600

plot_grid(A,B,nrow = 2,align = 'hv', rel_heights = c(1.5,2))
#1750x1050
