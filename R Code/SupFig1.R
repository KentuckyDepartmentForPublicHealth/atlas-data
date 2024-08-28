um<-read.csv("tsne_74_5000.csv",header = T)

plot(um[,c(2:3)], cex=0.5,xlab = "t-SNE 1", ylab = "t-SNE 2", main="FIt-SNE Plot (7375 samples)")

library(dbscan)

res <- optics(um[,c(2:3)], eps = 10, minPts = 3) #10
res <- extractXi(res, xi = 0.03) #0.025
hullplot(um[,c(2:3)],res, xlab = "t-SNE 1", ylab = "t-SNE 2",main="Convex Cluster Hulls (OPTICS)")
res$cluster

kNNdistplot(um[,c(2:3)], k = 2)
abline(h=1.75, col = "red", lty=2)
res <- dbscan(um[,c(2:3)], eps = 1.75, minPts = 3)
hullplot(um[,c(2:3)], res,xlab = "t-SNE 1", ylab = "t-SNE 2",main="Convex Cluster Hulls (DBSCAN)")
