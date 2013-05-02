kinome.grid<-function(X,title="")

{
	slide<-matrix(0,53,18)

	for (i in 1:18)
		{
		i1<-(i-1)*53+1
		i2<-i*53
		slide[,19-i]<-X$value[i1:i2]
		}
	slide[which(slide<1)]<-1

	n<-dim(slide)[1]
	p<-dim(slide)[2]

	maxi<-max(slide)/10

	plot(slide,xlim=range(1,p),ylim=range(1:n),type="n",
			xlab="columns",ylab="rows",main=title)	

	for (i in 1:p)
		for (j in 1:n)
			points(x=i,y=(54-j),pch=19,cex=log(slide[j,i]/maxi,base=2))

}

