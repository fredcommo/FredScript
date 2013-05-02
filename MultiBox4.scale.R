MultiBox4.scale<-function(class, val, xlab="class", ylab="values", nchar=10, Scale = TRUE, title=NULL, add.legend=TRUE, legends=NA, legcex=0.75)

# fonction graphique pour la représentation boxplot de plusieurs variables
# en fonction des niveaux de classe (1 facteur)

# A partir des données centrées réduites

{
	if (!is.factor(class)) class<-as.factor(class)
	
	val.scale <- as.data.frame(val)
	if (Scale) val.scale <- scale(as.data.frame(val))	

	p <- dim(val)[2]			# Nb variables
	p1 <- round(p/2,0)
	x <- p*length(levels(class))
	wex <- 1/x
	Outcex = 10/ncol(val.scale)
	Outpch = 19
	Maxi <- max(val.scale,na.rm=T)
	mini <- min(val.scale,na.rm=T)
	l <- length(levels(class))	# Nb niveaux dans la classe
	increm <- 0.01*(p-1)			# incrémentation: décalage des boîtes


	plot(val.scale[,1]~class, boxwex=wex, outcex = Outcex, outpch = Outpch, ylim=c(mini-1,Maxi+0.5), at=1:l-increm, axes=FALSE, xlab=xlab,
	ylab=ylab, names= substr(levels(class), start=1, stop=nchar), col=rainbow(p,1)[1], outcol=rainbow(p,1)[1], main=title)

	for (i in 2:p1){
		plot(val.scale[,i]~class, boxwex=wex, outcex = Outcex, outpch = Outpch, ylim=c(mini-1,Maxi+0.5), at=1:l-increm+wex*(i-1), add=TRUE, axes=FALSE, 
		col=rainbow(p,1)[i], outcol=rainbow(p,1)[i], xlab=NULL, ylab=NULL, names=substr(levels(class), start=1, stop=nchar))
		}

	plot(val.scale[,(p1+1)]~class, boxwex=wex, outcex = Outcex, outpch = Outpch, ylim=c(mini-1,Maxi+0.5), 
		at=1:l-increm+wex*p1, add=TRUE, col=rainbow(p,1)[(p1+1)], outcol=rainbow(p,1)[(p1+1)], xlab=NULL, ylab=NULL,
		names=substr(levels(class), start=1, stop=nchar))

	for (i in (p1+2):p){
		plot(val.scale[,i]~class, boxwex=wex, outcex = Outcex, outpch = Outpch, ylim=c(mini-1,Maxi+0.5), at=1:l-increm+wex*(i-1), add=TRUE, axes=FALSE,
		col=rainbow(p,1)[i], outcol=rainbow(p,1)[i], xlab=NULL, ylab=NULL, names=substr(levels(class), start=1, stop=nchar))
		}


	abline(h=0, col="blue", lty=3)

	if (add.legend){

		if (is.na(legends[1])) legends <- colnames(val.scale)
		else legends <- as.vector(legends)

		if (p<=20)	{	
			legend("bottomleft", legend=legends[1:p1], fill=rainbow(p,1)[1:p1], cex=legcex, bty="n")
			legend("bottomright", legend=legends[c(p1+1):p], fill=rainbow(p,1)[(p1+1):p], cex=legcex, bty="n")
			}

		else {
			p.leg <- floor(p/3)
			legend("bottomleft", legend=legends[1:p.leg], fill=rainbow(p,1)[1:p.leg], cex=legcex, bty="n")
			legend("bottom", legend=legends[(p.leg+1):(p.leg*2)], fill=rainbow(p,1)[(p.leg+1):(p.leg*2)], cex=legcex, bty="n")
			legend("bottomright", legend=legends[(p.leg*2+1):p], fill=rainbow(p,1)[(p.leg*2+1):p], cex=legcex, bty="n")
			}
		}
}




