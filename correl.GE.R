correl.GE<-function(Gene1, Gene2, Values){

# A corriger !!


	val<-Values

	gene1.Id<-names(which(Symb==Gene1))
	gene2.Id<-names(which(Symb==Gene2))

	gene1.values<-t(val[which(rownames(val) %in% gene1.Id),])
	gene2.values<-t(val[which(rownames(val) %in% gene2.Id),])

	ymin<-min(cbind(gene1.values, gene2.values))
	ymax<-max(cbind(gene1.values, gene2.values))

	C1<-hcl(240,seq(100,40,by=-10),seq(80,30,by=-5))	# 240~blue
	C2<-hcl(340,seq(100,40,by=-20),seq(60,10,by=-5))	# 350~red

	plot(gene2.values[,1]~gene1.values[,1], pch=19, col=C2[1], ylim=range(ymin,ymax))

	if (ncol(gene1.values)>1)
	for (i in 2:ncol(gene1.values))
		points(gene1.values[,i]~gene1.values[,1], pch=i, col=C1[i])

	if (ncol(gene2.values)>1)
	for (i in 2:ncol(gene2.values))	
		points(gene2.values[,i]~gene1.values[,1], pch=i, col=C2[i])

}

correl.GE("NOTCH4","NOTCH2",values)

