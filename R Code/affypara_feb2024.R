setwd("E:/data")
i<-read.csv("Atlas_Data.csv",header=T)

setwd("E:/data/AffyData")

library(affyPara)
library(hgu133plus2hsentrezgcdf)
library(hgu133plus2hsentrezg.db)

stopCluster(cl)
cl <-makeCluster(13, type="SOCK")

expdata <- preproPara(i$Filename, cdfname = "hgu133plus2hsentrezgcdf",
                      bgcorrect = TRUE, bgcorrect.method = "rma",
                      normalize = TRUE, normalize.method = "quantiles",
                      pmcorrect.method = "pmonly",
                      summary.method = "medianpolish",
                      verbose = T)

final <- exprs(expdata)
Annot <- data.frame(ENTREZID=sapply(contents(hgu133plus2hsentrezgENTREZID), paste, collapse=", "), chromosome=sapply(contents(hgu133plus2hsentrezgCHR), paste, collapse=", "))
all <- merge(Annot, final, by.x=0, by.y=0, all=T)
#Formatting dataframe, removing unnecessary column
rownames(all) <- all$Row.names
all$Row.names <- NULL

final<-NULL
expdata<-NULL

library(WGCNA)
maxmean_dat <- collapseRows(all[,3:dim(all)[2]], rowGroup = all$ENTREZID, rowID = rownames(all), method= "MaxMean", selectFewestMissing=FALSE)
entrezid <- row.names(maxmean_dat$datETcollapsed)
expression <- maxmean_dat$datETcollapsed
mydat.maxmean <- data.frame(entrezid, expression)
#mydat.maxmean <- head(mydat.maxmean,-1)
mydat.maxmean$entrezid<-NULL

maxmean_dat<-entrezid<-expression<-NULL

chromosomes<-all[,1:2]
mydat.maxmean.chr <- merge(chromosomes, mydat.maxmean, by.x=1, by.y=0, all=F)
mydat.maxmean.chr<-mydat.maxmean.chr[!duplicated(mydat.maxmean.chr), ]

Annot<-all<-chromosomes<-mydat.maxmean<-NULL

all<-mydat.maxmean.chr
rownames(all)<-all$ENTREZID
all$ENTREZID<-all$chromosome<-NULL

rownames(mydat.maxmean.chr)<-mydat.maxmean.chr$ENTREZID
mydat.maxmean.chr$ENTREZID<-mydat.maxmean.chr$chromosome<-NULL