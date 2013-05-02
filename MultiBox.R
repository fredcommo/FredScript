MultiBox<-function(class,val,xlab="class",ylab="values",nchar=4,title=NULL)

# fonction graphique pour la représentation boxplot de plusieurs variables
# en fonction des niveaux de classe
# class : vecteur class
# val: matrice de valeurs
# nchar : Nombre de caractères pris en compte pour le nom des boîtes (x.axis)

# A partir des données centrées réduites

{
	
	p<-dim(val)[2]			# Nb variables
	Maxi<-max(val,na.rm=T)
	mini<-min(val,na.rm=T)
	l<-length(levels(class))	# Nb niveaux dans la classe
	increm<-0.1*(p-1)			# incrémentation: décalage des boîtes


	plot(val[,1]~class,boxwex=0.15,ylim=c(mini-1,Maxi+0.5),at=1:l-increm,
			xlab=xlab,ylab=ylab,names= substr(levels(class),start=1,stop=nchar),
			col=rainbow(10,1)[1],main=title)

	for (i in 2:p)
	{
		plot(val[,i]~class,boxwex=0.15,ylim=c(mini-1,Maxi+0.5),at=1:l-increm+0.15*(i-1),add=TRUE,
			col=rainbow(10,1)[i],xlab=xlab,ylab=ylab,names=substr(levels(class),start=1,stop=nchar))
	}

	if (p<=3)	legend("topleft",colnames(val),fill=rainbow(10,1)[1:3],cex=0.75,bty="n")

	else	{	
		legend("topleft",colnames(val)[1:3],fill=rainbow(10,1)[1:3],cex=0.75,bty="n")
		legend("topright",colnames(val)[4:p],fill=rainbow(10,1)[4:p],cex=0.75,bty="n")
		}
}





