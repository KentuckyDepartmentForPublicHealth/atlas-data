library(smotefamily)
set.seed(1)

####Select identified samples

anno<-read.csv("Atlas_Data.csv",header = T)
anno<-anno[anno$Tumor!="#N/A",]

expdata<-readRDS("ExpData.rds")
expdata<-expdata[,anno$filename]
expdata<-data.frame(t(expdata))
expdata$Hist<-anno$Tumor

####Generate synthetic samples

maximum=966
x<-as.data.frame(summary(factor(anno$Tumor)))
x$total<-x$`summary(factor(anno$Tumor))`
x$`summary(factor(anno$Tumor))`<-NULL
x$tumor<-rownames(x)
x<-x[x$total<maximum,]
x<-x[order(x$total),]

syn = expdata[FALSE,]
for(i in 1:length(x$tumor)){
  s<-SMOTE(expdata[,-20361], expdata[,20361], K=min(x[i,1]-1,10), dup_size = ceiling(maximum/x[i,1]))
  expdata<-expdata[expdata$Hist!=x[i,2],]
  d<-s$syn_data
  d<-d[1:(maximum-x[i,1]),]
  syn<-rbind(syn,d)
}
colnames(syn)[20361]<-"Hist"
summary(factor(syn$Hist))

expdata<-readRDS("ExpData.rds")
expdata<-expdata[,anno$filename]
expdata<-data.frame(t(expdata))
expdata$Hist<-anno$Tumor

final<-data.frame(rbind(expdata,syn))

final<-readRDS("Final_SMOTE_dataset.RDS")

####final is balanced synthetic data with 966 samples/diagnosis; build a CV random forest model

library(ranger)

#grid search with ranger itself to select min.node.size and num.trees (the latter can't be tuned with caret with CV-fold)
hyper_grid <- expand.grid(
  mtry = c(170,175,180),
  min.node.size = c(1,2,3,4,5,6,7,8,9,10)
)

for(i in seq_len(nrow(hyper_grid))) {
  # fit model for ith hyperparameter combination
  fit <- ranger(
    x=final[,-20361], 
    y=factor(make.names(final$Hist)),
    num.trees       = 1000,
    min.node.size   = hyper_grid$min.node.size[i],
    mtry            = hyper_grid$mtry[i],
    verbose         = T,
    seed            = 1,
  )
  # export OOB error 
  hyper_grid$rmse[i] <- sqrt(fit$prediction.error)
}

# assess top 10 models
library(dplyr)
hyper_grid %>%
  arrange(rmse) %>%
  mutate(perc_gain = (0.01 - rmse) / 0.01 * 100) %>%
  head(5)
#lowest was 1000

#mtry min.node.size        rmse  perc_gain
#1  175             5 0.008866586 11.3341372
#2  180             8 0.008866586 11.3341372
#3  175             1 0.009913145  0.8685518
#4  200             1 0.009913145  0.8685518
#5  175             9 0.009913145  0.8685518
#6  180             2 0.009913145  0.8685518
#7  175             4 0.009913145  0.8685518
#8  180             4 0.009913145  0.8685518
#9  180             6 0.009913145  0.8685518
#10  142             1 0.010859306 -8.5930607

library(caret)
#cross-validation with grid with caret 
tgrid <- expand.grid(
  .mtry = 175, #around sqrt(20361) 175 was the best of 142,170,175,180,200
  .splitrule = "gini",
  .min.node.size = 5) #min.node.size 5 was the best from 1 to 9

rf <- caret::train(x=final[,-20361], y=factor(make.names(final$Hist)), 
                   method='ranger', 
                   metric='Accuracy',
                   tuneGrid=tgrid, 
                   trControl=trainControl(method='repeatedcv', number=5, repeats=3, verboseIter = T),
                   num.trees=1000, #best num.trees from ranger grid
                   importance="permutation")

saveRDS(rf,"RandomForest.RDS")
rf<-readRDS("RandomFOrest.RDS")

anno<-read.csv("tsne_74_5000.csv",header = T)
anno<-anno[anno$Tumor=="#N/A",]
newsamples<-readRDS("ExpData.rds")
newsamples<-newsamples[,anno$filename]
newsamples<-t(newsamples)

anno<-read.csv("tsne_74_5000.csv",header = T)
anno<-anno[anno$Tumor!="#N/A",]

colnames(newsamples)<-make.names(colnames(newsamples))

pred<-predict(rf, newsamples)
pred<-as.numeric(as.character(gsub("X","",pred)))
y<-cbind(rownames(newsamples),(pred))
y<-cbind(rownames(newsamples),as.character(unique(sort(factor(anno$Tumor)))[pred+1]))

write.csv(y, "rf.csv")


#confusion matrix
library(caret)
trainedsamples<-readRDS("ExpData.rds")
trainedsamples<-trainedsamples[,anno$filename]
trainedsamples<-t(trainedsamples)
colnames(trainedsamples)<-make.names(colnames(trainedsamples))
pred<-predict(rf, trainedsamples)
pred<-as.numeric(as.character(gsub("X","",pred)))
y<-cbind(rownames(trainedsamples),(pred))
y<-cbind(rownames(trainedsamples),as.character(unique(sort(factor(anno$Tumor)))[pred+1]))

o<-as.factor(anno$Tumor)
p<-as.factor(as.character(unique(sort(factor(anno$Tumor)))[pred+1]))
conf<-confusionMatrix(as.factor(anno$Tumor), as.factor(as.character(unique(sort(factor(anno$Tumor)))[pred+1])))
which(o!=p)

saveRDS(p,"rf_predict.RDS")
saveRDS(conf,"rf_conf.RDS")



