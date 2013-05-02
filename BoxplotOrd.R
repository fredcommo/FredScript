# sub-Boxplost post-ANOVA
t.repair <- as.data.frame(t(Repair.eset[cut.keep, ]))
t.repair <- scale(t.repair, scale=T)
sub1<-as.data.frame(t.repair[km2$cluster==1,])
sub2<-as.data.frame(t.repair[km2$cluster==2,])

col1= "royalblue3"
col2= "indianred3"

Smed<-apply(sub1, 2, median)
ord<-order(Smed, decreasing = TRUE)
boxplot(sub1[,ord], col=0, border=0, names= NULL, pars=par(cex.axis=0.5), horizontal=F, xlab="Expr. Z-score", main="Significant probes\nordered by cluster1 medians", ylim=range(-5,5))
#abline(v=seq(10,200,by=10), lty=3, col="mistyrose")

boxplot(sub1[,ord], col=col1, border=col1, names=NA, horizontal=F, add=T, notch = T, outpch = 19, outcex=0.25)
boxplot(sub2[,ord], col=col2, border=col2, names=NA, horizontal=F, add=T, notch = T, outpch = 19, outcex=0.25)

legend("topleft", legend=c("cluster1","cluster2"), fill=c(col1,col2), cex=1, bty="n")

