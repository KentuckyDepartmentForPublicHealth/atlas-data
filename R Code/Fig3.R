d<-read.csv("inconst.csv",header=T) #4 classifier columns with predicts where at least 1 is different

#unique inconsistency between any two classifiers
l <- list() 
for (x in 1:dim(d)[1]) {
  a<-data.frame()
  
  a[1:3,1]<-d[x,1]
  a[1:3,2]<-d[x,2:4]
  a[4:5,1]<-d[x,2]
  a[4:5,2]<-d[x,3:4]
  a[6,1]<-d[x,3]
  a[6,2]<-d[x,4]

  l[[x]]<-a
}

library(dplyr)
all<-bind_rows(l)
incon<-all[all[,1]!=all[,2],] #remove consistencies

#Generate a frame where the final number includes the sum of symmetrical relationships of X-Y and Y-X)
f<-incon %>% group_by_all() %>% summarise(COUNT = n())
f$dup<-duplicated(t(apply(f[,1:2], 1, sort)))
unique<-f[f$dup=="FALSE",]
dups<-f[f$dup=="TRUE",]
dupsflipped<-dups
dupsflipped$V1<-dups$tsne
dupsflipped$tsne<-dups$V1
f2<-rbind(unique,dupsflipped)
f3<-f2 %>% group_by(V1,tsne) %>% summarise(COUNT = sum(COUNT))

library(tidyr)
final<-f3 |> uncount(COUNT)

#dotplot of normalized inconsistencies
anno<-read.csv("Atlas_Data.csv", header=T, row.names = 1, stringsAsFactors = T)
anno<-anno[anno$Tumor!="#N/A",]
anno<-anno[anno$Tumor!="Nerve Ganglia",]

t<-(summary(factor(c(final$V1,final$tsne))))/sum(summary(factor(c(final$V1,final$tsne))))
z<-t/(summary(factor(anno$Tumor))/sum(summary(factor(anno$Tumor))))
z<-as.data.frame(z)
z$Diagnosis<-rownames(z)

library(ggpubr)
A<-ggdotchart(z,"Diagnosis", "z",add="segments", sorting = "descending", 
           ylab = "Inconsistent Predictions Involving the Diagnosis\nNormalized By the Number of Training Samples")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#1000x550


#Generate a chord diagram
# Transform input data in a adjacency matrix
adjacencyData <- with(final, table(V1, tsne)) #this is the same as f3, which can be entered into ChordDiagram

# Charge the circlize library
library(circlize)
chordDiagram(adjacencyData, transparency = 0.5)

#chord for only PA, GG, DIG, PXA - top 4 from above
pa<-final[final$V1=="Angiocentric" | final$tsne=="Angiocentric" | 
          final$V1=="DIG" | final$tsne=="DIG" | 
          final$V1=="PXA" | final$tsne=="PXA" |
          final$V1=="Ganglioglioma" | final$tsne=="Ganglioglioma",]

pa<-read.csv("pa.csv", header = T, stringsAsFactors = T) #the excel file from above rearranged to have the above diagnosis in one column

adjacencyData <- with(pa, table(tsne,V1))

chordDiagram(adjacencyData, transparency = 0.5,  annotationTrack = c("name", "grid"),
             annotationTrackHeight = c(0.03, 0.05), big.gap = 30)




chordDiagram(adjacencyData, transparency = 0.5, annotationTrack = "grid", big.gap = 30,
             preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(adjacencyData))))))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important

#1200x1200
