MultiBox5.scale<-function(X, val, class, xlab="class", ylab="values", Scale=TRUE, y.range=NA, title=NULL, subtitle=NULL, legcex=1, legpos=NA, col.palette=NULL)

# fonction graphique pour la représentation boxplot d'une variable
# en fonction des niveaux de plusieurs classes (n facteur)

# A partir des données centrées réduites

{
	if (!is.factor(class)) class<-as.factor(class)
	
	val.scale <- as.numeric(val)
	if (Scale) val.scale<-scale(as.numeric(val))	

	p <- nlevels(class)				# Nb de niveaux
	p1 <- round(p/2,0)
	x <- p*nlevels(class)
	wex <- 1/x
	mini <- min(val.scale,na.rm=T)
	Maxi <- max(val.scale,na.rm=T)
	l <- nlevels(as.factor(X))			# Nb niveaux dans la classe
	increm <- 0.01*(p-1)				# incrémentation: décalage des boîtes

	if (length(y.range)>1){
		mini=y.range[1]
		Maxi=y.range[2]
	}

	boxplot(val.scale[class==levels(class)[1]]~X[class==levels(class)[1]], boxwex=wex, ylim=c(mini-0.5,Maxi+0.5), at=1:l, axes=TRUE, xlab=xlab,
	ylab=ylab, names=unique(X), border=NA, main = title, sub=subtitle)

	if (is.null(col.palette)) col.palette <- colors()[seq(11, (11+p)*5, by=5)]

	for (i in 1:p){
		boxplot(val.scale[class==levels(class)[i]]~X[class==levels(class)[i]],
			boxwex=wex, ylim=c(mini-0.5, Maxi+0.5), add=T, at=1:l-increm+wex*(p1-i), axes=FALSE, xlab=xlab,
			ylab=ylab, names=unique(X), col=col.palette[i], pars=list(medlwd=1), outpch=8, outcex=0.75, outcol=col.palette[i])	#rainbow(p,1)
		}

	if (Scale) abline(h=0,col="blue",lty=3)

	legends <- levels(class)

	if (is.na(legpos)){
		pos <- "topleft"
		if (which.max(val.scale)<which.min(val.scale)) pos <- "topright"	
		}
		else pos <- legpos
		legend(pos, legend=legends, fill=col.palette[1:p], cex=legcex, bty="n")	#rainbow(p,1)[1:p]
}




