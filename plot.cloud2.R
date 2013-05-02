plot.cloud2<-function(X, V.lim = NULL, H.lim = NULL, Frontiers = F, Pch = 19, Cex = 0.5, Cols = "orangered")

# 2DPlot of a cloud with its marginal densities.

# X : matrix or data.frame with 2 columns

{
	def.par <- par(no.readonly = TRUE) # save default, for resetting... 

	x <-X[,1];	name1<-colnames(X)[1]
	y <-X[,2];	name2<-colnames(X)[2]
 
	source("D:\\Stats\\Doc R\\Scripts R\\density.limits.R")
	BIC.x<-mclustBIC(x)
	x.model<-mclustModel(x,BIC.x)

	BIC.y<-mclustBIC(y)
	y.model<-mclustModel(y,BIC.y)

	V.lim<- density.limits(x.model)
	H.lim<- density.limits(y.model)

	xdens <- density(x) 
	ydens <- density(y) 
	
	#layout.show(nf)
	nf <- layout(matrix(c(2,0,1,3), 2, 2, byrow = TRUE), c(3,1), c(1,3), TRUE) 
	layout.show(nf)	# à vérifier

	par(mar=c(4.2,4.2,1,1)) 
	plot(y~x, xlab = name1, ylab = name2, pch = Pch, cex = Cex, col = Cols)

	if (Frontiers){
		abline(h=H.lim, lty=3, col="blue"); abline(v=V.lim, lty=3, col="blue")

		# frontiers coord
		for (i in 1:length(V.lim))
			legend(x=V.lim[i],y=max(y)*0.95,x.intersp=0,y.intersp=-2,legend=round(V.lim[i],0),bty="n",text.col="blue",cex=0.9)
		for (i in 1:length(H.lim))
			legend(x=max(x)*0.90,y=H.lim[i],x.intersp=0,y.intersp=-2,legend=round(H.lim[i],0),bty="n",text.col="blue",cex=0.9)
		}

	par(mar=c(0,4.2,1,1)) 
	plot(xdens, axes=FALSE, col="darkblue", main=paste(name1,"density",sep=" "),ylab="",xlab="") 
	
	par(mar=c(4.2,0,1,1)) 
	plot(ydens$x~ydens$y, axes=FALSE, type="l", col="darkblue", main=paste(name2,"density",sep=" "),ylab="",xlab="") 

	par(def.par) 	# default : par(mar=c(5,4,2,2))
}
