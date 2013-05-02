MultiBox2<-function(class, val, xlab="class", ylab="values", nchar=10, title=NULL, legcex=0.75)

# fonction graphique pour la représentation boxplot de plusieurs variables
# en fonction des niveaux de classe
# class : vecteur class
# val: matrice de valeurs
# nchar : Nombre de caractères pris en compte pour le nom des boîtes (x.axis)

# A partir des données centrées réduites

{
	if (!is.factor(class)) class<-as.factor(class)

	p<-dim(val)[2]			# Nb variables
	p1<-round(p/2,0)
	x<-p*length(levels(class))
	wex<-1/x				# was set at 0.15
	Maxi<-max(val,na.rm=T)
	mini<-min(val,na.rm=T)
	l<-length(levels(class))	# Nb niveaux dans la classe
	increm<-0.01*(p-1)			# incrémentation: décalage des boîtes


	plot(val[,1]~class,boxwex=wex,ylim=c(mini-1,Maxi+0.5),at=1:l-increm,
			xlab=xlab,ylab=ylab,names= substr(levels(class),start=1,stop=nchar),
			col=rainbow(p,1)[1], axes=FALSE, main=title)

	for (i in 2:p1)
	{
		plot(val[,i]~class,boxwex=wex,ylim=c(mini-1,Maxi+0.5),at=1:l-increm+wex*(i-1),add=TRUE, axes=FALSE,
			col=rainbow(p,1)[i],xlab=xlab,ylab=ylab,names=substr(levels(class),start=1,stop=nchar))
	}

		plot(val[,(p1+1)]~class,boxwex=wex,ylim=c(mini-1,Maxi+0.5),at=1:l-increm+wex*p1,add=TRUE,
			col=rainbow(p,1)[c(p1+1)],xlab=xlab,ylab=ylab,names=substr(levels(class),start=1,stop=nchar))

	for (i in (p1+1):p)
	{
		plot(val[,i]~class,boxwex=wex,ylim=c(mini-1,Maxi+0.5),at=1:l-increm+wex*(i-1),add=TRUE, axes=FALSE,
			col=rainbow(p,1)[i],xlab=xlab,ylab=ylab,names=substr(levels(class),start=1,stop=nchar))
	}



	if (p<=3)	legend("topleft",colnames(val),fill=rainbow(p,1)[1:3],cex=legcex,bty="n")

	else	{
		legend("topleft",colnames(val)[1:p1],fill=rainbow(p,1)[1:p1],cex=legcex, bty="n")
		legend("topright",colnames(val)[(p1+1):p],fill=rainbow(p,1)[(p1+1):p],cex=legcex, bty="n")
		}
}





