setwd("~/Library/CloudStorage/OneDrive-Personal/Research/Tumor - Basic")
ExpData<-readRDS("ExpData.rds")

###Use t-sne based diagnosis for combat 
{
  
setwd("~/Library/CloudStorage/OneDrive-uoflhealth/Atlas")
anno<-read.csv("Atlas_data.csv",header = T)

d<-anno[,c(1,21,22)]
rownames(d)<-d$Filename
d$Filename<-NULL
d<-as.matrix(d)
library(parallelDist)
dis<-parDist(d)
dismatfull<-as.matrix(dis)

core<-anno[anno$Diagnosis!="#N/A",]
coremat<-dismatfull[rownames(dismatfull) %in% core$Filename,colnames(dismatfull) %in% core$Filename]

mattest = matrix(, nrow = dim(coremat)[1], ncol = 4)
colnames(mattest)<- c("Filename", "closest","predict","orig")
mattest[,1]<-rownames(coremat) #add sample samples to test
for (i in 1:dim(coremat)[1]){
  mattest[i,2]<-colnames(coremat)[which(coremat[i,]==min(coremat[i,][coremat[i,]>0]))][1] #the 1 at the end is incase of more than one samples with the same distance
  mattest[i,3]<-anno[which(anno$Filename==mattest[i,2]),17]
  mattest[i,4]<-anno[which(anno$Filename==mattest[i,1]),17]
}
mattest<-data.frame(mattest, stringsAsFactors = T)

library(caret)
confusionMatrix(data=mattest$predict,reference = mattest$orig) #predicted vs. true statistics

test<-anno[anno$Diagnosis=="#N/A",]
testmat<-dismatfull[rownames(dismatfull) %in% test$Filename,!colnames(dismatfull) %in% test$Filename]

mattest = matrix(, nrow = dim(testmat)[1], ncol = 3)
colnames(mattest)<- c("Filename", "closest","predict")
mattest[,1]<-rownames(testmat) #add sample samples to test
for (i in 1:dim(testmat)[1]){
  mattest[i,2]<-colnames(testmat)[which(testmat[i,]==min(testmat[i,][testmat[i,]>0]))][1] #the 1 at the end is incase of more than one samples with the same distance
  mattest[i,3]<-anno[which(anno$Filename==mattest[i,2]),17]
}

}

diag<-read.csv("/Users/amistry/Library/CloudStorage/OneDrive-uoflhealth/Atlas/Atlas_data.csv")
diag<-diag[,c(1,3,4,6,7,23)]
rownames(diag)<-diag$Filename
diag$Filename<-NULL
colnames(diag)<-c("Year","Dataset","Institution","Country","Diagnosis")
meta<-diag[colnames(ExpData),]

#PCVA
library("ExpressionNormalizationWorkflow")
inpData <- expSetobj(ExpData, meta)
cvrts_eff_var <- inpData@phenoData@varMetadata[["labelDescription"]]
rm(diag,meta,ExpData)
pvcAnaly(inpData, 0.5, cvrts_eff_var)
#pvcAnaly(inpData, 0.75, cvrts_eff_var)


library(ggpubr)
library(scales)
prey <- c(0.0904312, 0.1644299, 0.4682778, 0.01516156, 0.05220623, 0.2094933)
x_labels <- c("Dataset", "Institution", "Diagnosis", "Country", "Year", "Residual")
data <- data.frame(Variable = x_labels, Percentage = prey)
data$Variable <- factor(data$Variable, levels = data$Variable[order(-data$Percentage)])
ggplot(data, aes(x = Variable, y = Percentage)) +
  geom_bar(stat = "identity") +
  labs(x = "Variable", y = "Percent of the Variance Explained", title = "PCVA") +
  scale_y_continuous(labels = scales::percent) + theme_pubr() +
  geom_text(aes(label = scales::percent(Percentage)), vjust = -0.5, color = "black", size = 5)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

#combat
library(sva)
library(BiocParallel)

#adjust sequentially all
serialParam <- SerialParam()
modcombat = model.matrix(~as.factor(Diagnosis),data=meta)

batch<-as.factor(meta$Institution)
combat_edata = ComBat(dat=ExpData, batch=batch, mod=modcombat, 
                      par.prior=T, prior.plot=F,ref.batch = NULL,
                      BPPARAM = bpparam("SerialParam"))

batch<-as.factor(meta$Dataset)
combat_edata = ComBat(dat=combat_edata, batch=batch, mod=modcombat, 
                      par.prior=T, prior.plot=F,ref.batch = NULL,
                      BPPARAM = bpparam("SerialParam"))

batch<-as.factor(meta$Year)
combat_edata = ComBat(dat=combat_edata, batch=batch, mod=modcombat, 
                      par.prior=T, prior.plot=F,ref.batch = NULL,
                      BPPARAM = bpparam("SerialParam"))

NewExp<-as.data.frame(combat_edata)
#saveRDS(NewExp, "ExpData_BC.RDS")

set.seed(1)
t<-fftRtsne(t(combat_edata))

combat<-as.data.frame(colnames(combat_edata))
combat<-cbind(combat, t)
rownames(combat)<-combat$`colnames(combat_edata)`
combat$`colnames(combat_edata)`<-NULL

d<-read.csv("/Users/amistry/Library/CloudStorage/OneDrive-uoflhealth/Atlas/Atlas_data.csv", header=T, row.names = 1, stringsAsFactors = T)
new<-merge(combat,d,by="row.names",all.x=TRUE)
rownames(new)<-new$Row.names
new$Row.names<-NULL
new$V1<-new$V2<-NULL

colnames(new)[1]<-"V1"
colnames(new)[2]<-"V2"

d<-new
d<-d[,c(1,2,24)]
d$Final<-d$combat_diagnosis

library(ggpubr)
library(plotly)
p<-ggplot()+geom_point(data=d, aes(V1,V2,color=factor(Final), label=rownames(d)), size=0.75)+theme_pubr()+
  theme(legend.position = "right")+rremove("axis")+rremove("axis.text")+rremove("ticks")+rremove("xylab")+
  theme(panel.border = element_rect(colour = "black", fill=NA))
ggplotly(p)

library("ExpressionNormalizationWorkflow")
inpData <- expSetobj(combat_edata, meta)
cvrts_eff_var <- inpData@phenoData@varMetadata[["labelDescription"]]
rm(diag,meta,ExpData)
pvcAnaly(inpData, 0.5, cvrts_eff_var)

library(ggpubr)
library(scales)
posty <- c(prey,0.01624525, 0.02577943, 0.6672807, 0.003763693, 0.004703307, 0.2822276)
postx <- rep(c("Dataset","Institution", "Temp. Diagnosis","Country","Year","Residual"),2)
Adjustment<-c(rep("Before",6),rep("After",6))
data <- data.frame(Variable = postx, Percentage = posty, Adjustment=Adjustment)
data$Variable <- factor(data$Variable, levels = c("Temp. Diagnosis","Institution","Dataset","Year","Country","Residual"))
data$Adjustment <- factor(data$Adjustment, levels = c("Before","After"))
ggplot(data, aes(x = Variable, y = Percentage, fill = Adjustment)) +
  geom_col(position = "dodge") +  # Prevent stacking by using "dodge" position
  labs(x = "Variable", y = "Percent of the Variance Explained", title = "Principal Variance Component Analysis", fill = "SVA Adjustment") +
  scale_y_continuous(labels = scales::percent) + 
  theme_pubr() +
  geom_text(aes(label = scales::percent(Percentage, accuracy = 0.1)),
            position = position_dodge(width = 0.9), vjust = -0.3, color = "black", size = 5) 
#900 x 500