kinome.grid2<-function(X,title="",coef=2000,log.base=10)

{
	slide<-matrix(0,53,18)

	for (i in 1:18){
		i1<-(i-1)*53+1
		i2<-i*53
		slide[,19-i]<-X$value[i1:i2]
		}
	

	n<-dim(slide)[1]
	p<-dim(slide)[2]

	plot(slide,xlim=range(1,p),ylim=range(1:n),type="n",
			xlab="columns",ylab="rows",main=title)	
	
	#version4 (voir programme R)
	slide[which(slide<0)]<-1

	for (i in 1:p)
		for (j in 1:n)
		
		#version4
		{	
		val<-ceiling(slide[j,i]/coef)
		points(x=i,y=(54-j),pch=19,cex=log(val,base=log.base))
		}
	
	# Repères
	for (j in 1:53){
		points(x=1,y=j,pch=19,cex=1,col="lightblue")
		points(x=18,y=j,pch=19,cex=1,col="lightblue")
	}
	for (j in c(1,14,27,40,53))
		for (i in seq(1,18))
			points(x=i,y=j,pch=19,cex=1,col="lightblue")
}
