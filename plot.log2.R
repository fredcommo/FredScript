plot.log2<-function(X,Chr.num)

	# plot.log2(aCGH.obj,Chr.num) représente les valeurs log2 de chaque patient pour un Chr passé en argument
	# positionne le centromere selon les données de human.chrom.info.Jul03
	# ajoute la moyenne par sonde

	# aCGH.obj : objet de classe aCGH
	# Chr.num : n° du chromosome à représenter

{
	log2<-log2.ratios(X)
	index<-which(clones.info(X)$Chrom==Chr.num)
	log2.Chr<-log2[index,]

	clones.kb<-clones.info(X)$kb[index]

	max<-max(log2.Chr)
	min<-min(log2.Chr)
	mean<-apply(log2.Chr,1,mean)

	centro<-human.chrom.info.Jul03$centromere[Chr.num]
	plot(log2.Chr[,1]~clones.kb,pch=20,cex=0.2,col="lightblue",ylim=range(min,max),
	ylab="Log2",xlab="probes pos kb",main=paste("Chrom",Chr.num,sep=" "))

	for (i in 2:dim(log2.Chr)[2])
		{
			points(log2.Chr[,i]~clones.kb,pch=20,cex=0.2,col="lightblue")
		}
	points(mean~clones.kb,pch=8,col="blue",cex=0.2)
	abline(v=centro,lty=3,col="red")
	legend("topleft",legend="mean",bty="n",lty=1,col="blue",cex=0.75)
}
	