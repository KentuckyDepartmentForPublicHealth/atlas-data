library(smotefamily)
library(plotly)
library(ggplot2)
set.seed(1)

anno<-read.csv("/Users/amistry/Desktop/tsne_74_5000.csv",header = T)

expdata<-readRDS("/Users/amistry/Library/CloudStorage/OneDrive-Personal/Research/Tumor - Basic/ExpData.rds")
expdata<-data.frame(t(expdata[,anno$filename]))
expdata$Hist<-anno$classifier
expdata<-expdata[expdata$Hist!="Met - Colon",]
expdata<-expdata[expdata$Hist!="Met - Esophageal",]

barplot(summary(reorder(factor(expdata$Hist),factor(expdata$Hist), FUN=length)), las=2)

maximum=966
x<-as.data.frame(summary(factor(anno$classifier)))
x$total<-x$`summary(factor(anno$classifier))`
x$`summary(factor(anno$classifier))`<-NULL
x$tumor<-rownames(x)
x<-x[x$total<maximum,]
x<-x[order(x$total),]
x<-x[-c(1:2),]

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

expdata<-readRDS("/Users/amistry/Library/CloudStorage/OneDrive-Personal/Research/Tumor - Basic/ExpData.rds")
expdata<-data.frame(t(expdata[,anno$filename]))
expdata$Hist<-anno$classifier

final<-rbind(expdata,syn)

#see on t-SNE where the fake samples are
#syn$Hist <- sub("^", "SYN-", syn$Hist)
#setwd()<-"/Users/amistry/Desktop/FIt-SNE-master"
#source("fast_tsne.R")
#ts<-fftRtsne(as.matrix(final[,-20361]))
#tsne<-data.frame(ts)
#tsne$final<-final[,20361]
#ggplotly(ggplot(tsne[tsne$final!="#N/A",], aes(X1,X2,color=final))+
#           geom_point(size=1)+theme_classic())
#final$Hist <- sub("SYN-", "", final$Hist)

summary(reorder(factor(final$Hist),factor(final$Hist), FUN=length))
final<-final[final$Hist!="Met - Colon",]
final<-final[final$Hist!="Met - Esophageal",]

library(dplyr)
traintest<-final[final$Hist!="#N/A",]
traintest$Hist<-factor(traintest$Hist)
traintest$filename<-rownames(traintest)
test<-traintest[traintest$Hist!="#N/A",] %>% group_by(Hist) %>% slice_sample(prop=0.2)
train<-traintest[!traintest$filename %in% test$filename,]
test<-traintest[!traintest$filename %in% train$filename,]
test$filename<-train$filename<-traintest$filename<-NULL
classify<-final[final$Hist=="#N/A",]

alltest<-final[final$Hist!="#N/A",]
s<-summary(reorder(factor(alltest$Hist),factor(alltest$Hist), FUN=length))
wts<-nrow(alltest)/(length(unique(alltest$Hist))*s)

anno<-d<-expdata<-s<-syn<-traintest<-x<-NULL

library(randomForest)
library(ranger)
library(caret)
# This fails for large features
#rf <- randomForest(Hist ~ ., data=data)
# This works
rf<-ranger(x=alltest[,-20361], y=factor(alltest[,20361]), importance = "permutation", scale.permutation.importance =T,
           oob.error = T, classification = T, num.trees = 2000, probability = T)
print(rf)


tgrid <- expand.grid(
  .mtry = c(142, 175, 200, 250),
  .splitrule = "gini",
  .min.node.size = 1)

library(ROSE)
library(themis)
rf <- caret::train(x=alltest[,-20361], y=factor(make.names(alltest$Hist)), 
            method='ranger', 
            metric='Accuracy',
            tuneGrid=tgrid, 
            trControl=trainControl(method='repeatedcv', number=5, repeats=3, verboseIter = T, 
                                   classProbs = T),
            num.trees=1000,
            importance="permutation")

p3<-predict(rf, classify[,-20361])

default_rmse <- sqrt(rf$prediction.error)
n_features<-20360

hyper_grid <- expand.grid(
  mtry = c(75,100,125,150,175,200),
  min.node.size = c(1, 3, 5, 10),
  replace = c(TRUE, FALSE),                               
  sample.fraction = c(.5, .63, .8),
  num.trees = c(250,500,1000),
  rmse = NA                                               
)

for(i in seq_len(nrow(hyper_grid))) {
  # fit model for ith hyperparameter combination
  fit <- ranger(
    x=train[,-20361], 
    y=train[,20361],
    num.trees       = hyper_grid$num.trees[i],
    mtry            = hyper_grid$mtry[i],
    min.node.size   = hyper_grid$min.node.size[i],
    replace         = hyper_grid$replace[i],
    sample.fraction = hyper_grid$sample.fraction[i],
    verbose         = T,
    seed            = 1
  )
  # export OOB error 
  hyper_grid$rmse[i] <- sqrt(fit$prediction.error)
}

hyper_grid %>%
  arrange(rmse) %>%
  mutate(perc_gain = (default_rmse - rmse) / default_rmse * 100) %>%
  head(10)

p1<-predict(rf,alltest[,-20361])
confusionMatrix(p1$predictions,alltest$Hist)

p2<-predict(rf,test[,-20361])
confusionMatrix(p2$predictions,test$Hist)

p3<-predict(rf,classify[,-20361], type = "response")
confusionMatrix(p3$predictions,classify$Hist)

classify$pred<-p3$predictions

p<-p3$predictions
rownames(p)<-rownames(classify)

mat = matrix(, nrow = dim(p)[1], ncol = 2)
colnames(mat)<- c("test", "predict")
mat[,1]<-rownames(p) #add sample samples to test
for (i in 1:dim(p)[1]){
  mat[i,2]<-colnames(p)[which(p[i,]==max(p[i,][p[i,]>0]))][1]
}

newpreds<-classify[,20361:20362]
write.csv(mat, "newpreds.csv")
