anno<-read.csv("tsne_74_5000.csv",header = T)

#####generate distance matrix for 'core' dataset
#library(parallelDist)

#anno<-anno[anno$Tumor!="#N/A",]

#expdata<-readRDS("/Users/amistry/Library/CloudStorage/OneDrive-Personal/Research/Tumor - Basic/ExpData.rds")
#expdata<-expdata[,anno$filename]
#expdata<-data.frame(t(expdata))

#d<-as.matrix(expdata)
#dis<-parDist(d)
#dismat<-as.matrix(dis)
#saveRDS(dismat, "dismat_core.RDS")
#dismat<-readRDS("dismat_core.RDS")

#generate distance matrix for all data ('core' and test)
#expdata<-readRDS("/Users/amistry/Library/CloudStorage/OneDrive-Personal/Research/Tumor - Basic/ExpData.rds")
#expdata<-expdata[,anno$filename]
#expdata<-data.frame(t(expdata))
#d<-as.matrix(expdata)
#dis<-parDist(d)
#dismatfull<-as.matrix(dis)

dismat<-readRDS("dismat_core.RDS") #'core' dataset
dismatfull<-readRDS("dismat_full.RDS") #includes 'core' and 'test' datasets

#max of the min distances
mins<-c()
for (i in 1:dim(dismat)[1]){
  mins[i]<-min(dismat[i,][dismat[i,]>0])
}
max(mins) #maximum "minimum" (non-zero) distance 

#closest match
mat = matrix(, nrow = dim(dismat)[1], ncol = 6)
colnames(mat)<- c("test", "closest", "distance", "close?","predict","orig")
mat[,1]<-rownames(dismat) #add sample samples to test
for (i in 1:dim(dismat)[1]){
  mat[i,2]<-colnames(dismat)[which(dismat[i,]==min(dismat[i,][dismat[i,]>0]))][1] #the 1 at the end is incase of more than one samples with the same distance
  mat[i,3]<-min(dismat[i,][dismat[i,]>0])
  mat[i,4]<-min(dismat[i,][dismat[i,]>0])<=max(mins)
  mat[i,5]<-anno[which(anno$filename==mat[i,2]),16]
  mat[i,6]<-anno[which(anno$filename==mat[i,1]),16]
}
mat<-data.frame(mat, stringsAsFactors = T)

library(caret)
dclass<-confusionMatrix(data=mat$predict,reference = mat$orig) #predicted vs. true statistics

library(ggpubr)
library(ggplot2)
#confusion matrix
prettyConfused<-function(Actual,Predict,colors=c("white","blue","red"),text.scl=1){
  actual = as.data.frame(table(Actual))
  names(actual) = c("Actual","ActualFreq")
  #build confusion matrix
  confusion = as.data.frame(table(Actual, Predict))
  names(confusion) = c("Actual","Predicted","Freq")
  #calculate percentage of test cases based on actual frequency
  confusion = merge(confusion, actual, by=c('Actual','Actual'))
  confusion$Percent = confusion$Freq/confusion$ActualFreq*100
  confusion$ColorScale<-confusion$Percent*-1
  confusion[which(confusion$Actual==confusion$Predicted),]$ColorScale<-confusion[which(confusion$Actual==confusion$Predicted),]$ColorScale*-1
  confusion$Label<-paste(round(confusion$Percent,0),"% (",confusion$Freq,")",sep="")
  confusionhighlight<-confusion[confusion$Percent!=0,]
  confusionhighlight<-confusionhighlight[confusionhighlight$Actual!=confusionhighlight$Predicted,]
  tile <- ggplot() +
    geom_tile(aes(x=Actual, y=Predicted,fill=ColorScale),data=confusion, color="lightgrey",size=0.1) +
    #geom_tile(aes(x=Actual, y=Predicted,fill=ColorScale),data=confusionhighlight, color="black",size=0.2) + #to border the cells 
    labs(x="Actual",y="Predicted")
  tile = tile +
    geom_text(aes(x=Actual,y=Predicted, label=Freq),data=confusionhighlight, size=3, colour="black") +
    scale_fill_gradient2(low=colors[2],high=colors[3],mid=colors[1],midpoint = 0, name="Percent")+ #add guide="none" to remove legend
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  tile
}

library(grid)
#x_grob <- rectGrob(x=1:52, y=0, gp=gpar(color='black', fill=rainbow(52), alpha=0.75))
#y_grob <- rectGrob(y=1:52, x=0, gp=gpar(color='black', fill=rainbow(52), alpha=0.75))
e<-prettyConfused(mat$orig, mat$predict)+coord_cartesian(clip='off')+
  annotate("text", x=1, y=52, label= "Euclidean Classifier: Overall Accuracy: 99.58% [95% CI: 99.37% to 99.74%]", hjust = 0, fontface = "bold")
#  annotation_custom(grob=x_grob, xmin=0, xmax=1, ymin=-3.75, ymax=3.75)+
#  annotation_custom(grob=y_grob, xmin=-3.75, xmax=3.75, ymin=0, ymax=1)
#750x750

####predict test samples
testmat<-dismatfull[!rownames(dismatfull) %in% rownames(dismat),colnames(dismatfull) %in% colnames(dismat)]

mattest = matrix(, nrow = dim(testmat)[1], ncol = 4)
colnames(mattest)<- c("filename", "closest","meets threshold?","predict")
mattest[,1]<-rownames(testmat) #add sample samples to test
for (i in 1:dim(testmat)[1]){
  mattest[i,2]<-colnames(testmat)[which(testmat[i,]==min(testmat[i,][testmat[i,]>0]))][1] #the 1 at the end is incase of more than one samples with the same distance
  mattest[i,3]<-min(testmat[i,][testmat[i,]>0])<=max(mins)
  mattest[i,4]<-anno[which(anno$filename==mattest[i,2]),16]
  }

#write.csv(mattest,"newdist.csv")

####t-sne based dismat####
anno<-read.csv("tsne_74_5000.csv",header = T)

d<-anno[,c(1:3)]
rownames(d)<-d$filename
d$filename<-NULL
d<-as.matrix(d)
library(parallelDist)
dis<-parDist(d)
dismatfull<-as.matrix(dis)

core<-anno[anno$Tumor!="#N/A",]
coremat<-dismatfull[rownames(dismatfull) %in% core$filename,colnames(dismatfull) %in% core$filename]

mattest = matrix(, nrow = dim(coremat)[1], ncol = 4)
colnames(mattest)<- c("filename", "closest","predict","orig")
mattest[,1]<-rownames(coremat) #add sample samples to test
for (i in 1:dim(coremat)[1]){
  mattest[i,2]<-colnames(coremat)[which(coremat[i,]==min(coremat[i,][coremat[i,]>0]))][1] #the 1 at the end is incase of more than one samples with the same distance
  mattest[i,3]<-anno[which(anno$filename==mattest[i,2]),16]
  mattest[i,4]<-anno[which(anno$filename==mattest[i,1]),16]
}
mattest<-data.frame(mattest, stringsAsFactors = T)

library(caret)
tclass<-confusionMatrix(data=mattest$predict,reference = mattest$orig) #predicted vs. true statistics

library(grid)
#x_grob <- rectGrob(x=1:52, y=0, gp=gpar(color='black', fill=rainbow(52), alpha=0.75))
#y_grob <- rectGrob(y=1:52, x=0, gp=gpar(color='black', fill=rainbow(52), alpha=0.75))
t<-prettyConfused(mattest$orig, mattest$predict)+coord_cartesian(clip='off')+
  annotate("text", x=1, y=52, label= "t-SNE-Based Classifier: Overall Accuracy: 99.92% [95% CI: 99.80% to 99.98%]", hjust = 0, fontface = "bold")
#  annotation_custom(grob=x_grob, xmin=0, xmax=1, ymin=-3.75, ymax=3.75)+
#  annotation_custom(grob=y_grob, xmin=-3.75, xmax=3.75, ymin=0, ymax=1)
#750x750

test<-anno[anno$Tumor=="#N/A",]
testmat<-dismatfull[rownames(dismatfull) %in% test$filename,!colnames(dismatfull) %in% test$filename]

mattest = matrix(, nrow = dim(testmat)[1], ncol = 3)
colnames(mattest)<- c("filename", "closest","predict")
mattest[,1]<-rownames(testmat) #add sample samples to test
for (i in 1:dim(testmat)[1]){
  mattest[i,2]<-colnames(testmat)[which(testmat[i,]==min(testmat[i,][testmat[i,]>0]))][1] #the 1 at the end is incase of more than one samples with the same distance
  mattest[i,3]<-anno[which(anno$filename==mattest[i,2]),16]
}

#write.csv(mattest,"newtsne.csv")

#GBM
pred<-readRDS("gbm_predict.RDS")
pred<-gsub("RB","Retinoblastoma",pred)
pred<-gsub("NB","Neuroblastoma",pred)
pred<-gsub("Neuroblastoma-FOXR2","NB-FOXR2",pred)
pred<-gsub("Pit. Tumor","PitNET",pred)

anno<-read.csv("tsne_74_5000.csv",header = T)
anno<-anno[anno$Tumor!="#N/A",]
g<-prettyConfused(as.factor(anno$Tumor), pred)+coord_cartesian(clip='off')+
  annotate("text", x=1, y=52, label= "LightGBM Classifier: Overall Accuracy: 99.98% [95% CI: 99.89% to 100%]", hjust = 0, fontface = "bold")

gbmconf<-readRDS("gbm_conf.RDS")


#RF
pred<-readRDS("rf_predict.RDS")
pred<-gsub("RB","Retinoblastoma",pred)
pred<-gsub("NB","Neuroblastoma",pred)
pred<-gsub("Neuroblastoma-FOXR2","NB-FOXR2",pred)
pred<-gsub("Pit. Tumor","PitNET",pred)

anno<-read.csv("tsne_74_5000.csv",header = T)
anno<-anno[anno$Tumor!="#N/A",]
r<-prettyConfused(as.factor(anno$Tumor), pred)+coord_cartesian(clip='off')+
  annotate("text", x=1, y=52, label= "Random Forest Classifier: Overall Accuracy: 100% [95% CI: 99.93% to 100%]", hjust = 0, fontface = "bold")

rfconf<-readRDS("rf_conf.RDS")


library(cowplot)
plot_grid(r,g,e,t, align = c("hv"), ncol=2, labels = 'AUTO', label_size = 20)
#1700x1500


ss = matrix(, nrow = 52*4, ncol = 4)
colnames(ss)<-c("Sensitivity","Specificity","Class","Classifier")
ss[1:52,1:2]<-dclass[["byClass"]][,1:2]
ss[1:52,3]<-rownames(dclass[["byClass"]])
ss[1:52,4]<-"Euclidean"

ss[53:104,1:2]<-tclass[["byClass"]][,1:2]
ss[53:104,3]<-rownames(tclass[["byClass"]])
ss[53:104,4]<-"t-SNE-Based"

ss[105:156,1:2]<-gbmconf[["byClass"]][,1:2]
ss[105:156,3]<-rownames(gbmconf[["byClass"]])
ss[105:156,4]<-"LightGBM"
  
ss<-data.frame(ss, stringsAsFactors = T)
ss$Class<-gsub('Class: ','',ss$Class)
ss$Class<-gsub("RB","Retinoblastoma",ss$Class)
ss$Class<-gsub("NB","Neuroblastoma",ss$Class)
ss$Class<-gsub("Neuroblastoma-FOXR2","NB-FOXR2",ss$Class)
ss$Class<-gsub("Pit. Tumor","PitNET",ss$Class)

ss$Class<-as.factor(ss$Class)
ss$Sensitivity<-as.numeric(as.character(ss$Sensitivity))
ss$Specificity<-as.numeric(as.character(ss$Specificity))

library(ggpubr)
ggdotchart(ss,"Class","Sensitivity",color="Classifier",shape="Classifier",add="segments")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggdotchart(ss,"Class","Specificity",color="Classifier",shape="Classifier",add="segments")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
