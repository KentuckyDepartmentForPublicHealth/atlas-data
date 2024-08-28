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

####final is balanced synthetic data with 966 samples/diagnosis; build a CV random forest model
library(caret)
library(lightgbm)
#split data
set.seed(1)
indexes = createDataPartition(final$Hist, p = .85, list = F)
final$Hist<-as.numeric(as.factor(final$Hist))-1

train = as.matrix(final[indexes, ])
test = as.matrix(final[-indexes, ])

train_x = train[, -20361]
train_y = train[, 20361]

test_x = test[, -20361]
test_y = test[, 20361]

dtrain = lgb.Dataset(train_x, label = train_y, params = list(max_bin=255), free_raw_data = FALSE) #change to 375 and 500
dtest = lgb.Dataset.create.valid(dtrain, data = test_x, label = test_y)

# validation data
valids = list(test = dtest)

grid_search <- expand.grid(boost = c("dart","gbdt"),
                           iter = c(100),
                           leaves = c(30,40,50,60),
                           rate =c(0.1,0.05,0.01))
model <- list()
perf <- numeric(nrow(grid_search))

for (i in 1:nrow(grid_search)) {
  model[[i]] <- lgb.train(list(objective="multiclass",
                               metric="multi_error",
                               num_class=52,   #52 if "cerebellum,fetal" added. 
                               seed=1,
                               boosting=grid_search[i, "boost"],
                               num_iterations=grid_search[i, "iter"],
                               num_leaves=grid_search[i, "leaves"],
                               learning_rate=grid_search[i, "rate"]),
                          dtrain,
                          nrounds=100,  #change this for the final model to 100
                          early_stopping_rounds = 10,
                          valids)
  perf[i] <- min(as.numeric(model[[i]]$record_evals$test$multi_error$eval)) #change this
}
cat("Model ", which.min(perf), " is multi_error: ", min(perf), sep = "","\n")
print(grid_search[which.min(perf), ])

#July --max_bin 255: best error 0.0002636436 with gbdt; 30 leaves; 0.1 learning rate
#July --max_bin 375: best error 0.0002636436 with gbdt; 30 leaves; 0.1 learning rate
#July --max_bin 500: best error 0.0003954653 with gbdt; 50 leaves; 0.05 learning rate

#Bin 255 - best error 0.000268745
#-best leaves 50 and rate 0.1 (best iteration is 63 for dart and 60 gbdt); max_bin 375, learning rate 0.2, or leaves 80 did not improve error rate

#building a cross-validated model based on above parameters
dtrain = lgb.Dataset(train_x, label = train_y, params = list(max_bin=255), free_raw_data = FALSE) #change to 375 and 500
dtest = lgb.Dataset.create.valid(dtrain, data = test_x, label = test_y)
valids = list(test = dtest)

finalmodel = lgb.train(list(objective="multiclass",
                       metric="multi_error",
                       num_class=52, #52 if "cerebellum,fetal" added. 
                       seed=1,
                       boosting="gbdt",
                       num_iterations=100,
                       num_leaves=30,
                       learning_rate=0.1),
                  dtrain,
                  nrounds=100,
                  early_stopping_rounds = 10,
                  valids)

anno<-read.csv("tsne_74_5000.csv",header = T)
anno<-anno[anno$Tumor=="#N/A",]
newsamples<-readRDS("ExpData.rds")
newsamples<-newsamples[,anno$filename]
newsamples<-t(newsamples)

anno<-read.csv("tsne_74_5000.csv",header = T)
anno<-anno[anno$Tumor!="#N/A",]

pred<-predict(finalmodel, newsamples)
pred_y = max.col(pred)-1
y<-cbind(rownames(pred),as.character(unique(sort(factor(anno$Tumor)))[pred_y+1]))
write.csv(y,"gbm.csv")

#confusion matrix
trainedsamples<-readRDS("ExpData.rds")
trainedsamples<-trainedsamples[,anno$filename]
trainedsamples<-t(trainedsamples)
pred<-predict(finalmodel, trainedsamples)
pred_y = max.col(pred)-1

o<-as.factor(anno$Tumor)
p<-as.factor(as.character(unique(sort(factor(anno$Tumor)))[pred_y+1]))
confusionMatrix(as.factor(anno$Tumor), as.factor(as.character(unique(sort(factor(anno$Tumor)))[pred_y+1])))
which(o!=p)

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
  confusion$Label<-paste(round(confusion$Percent,0),"%, n=",confusion$Freq,sep="")
  tile <- ggplot() +
    geom_tile(aes(x=Actual, y=Predicted,fill=ColorScale),data=confusion, color="lightgrey",size=0.1) +
    labs(x="Actual",y="Predicted")
  tile = tile +
    geom_text(aes(x=Actual,y=Predicted, label=""),data=confusion, size=text.scl, colour="black") +
    scale_fill_gradient2(low=colors[2],high=colors[3],mid=colors[1],midpoint = 0,guide='none')+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  tile
}

library(grid)
x_grob <- rectGrob(x=1:51, y=0, gp=gpar(color='black', fill=rainbow(51), alpha=0.75))
y_grob <- rectGrob(y=1:51, x=0, gp=gpar(color='black', fill=rainbow(51), alpha=0.75))
prettyConfused(as.factor(anno$Tumor), as.factor(as.character(unique(sort(factor(anno$Tumor)))[pred_y+1])))+coord_cartesian(clip='off')+
  annotation_custom(grob=x_grob, xmin=0, xmax=1, ymin=-3.75, ymax=3.75)+
  annotation_custom(grob=y_grob, xmin=-3.75, xmax=3.75, ymin=0, ymax=1)
#750x750