MultiBox.scale<-function(class,val,xlab="class",ylab="values",nchar=4,title=NULL)

# fonction graphique pour la représentation boxplot de plusieurs variables
# en fonction des niveaux de classe

# A partir des données centrées réduites

{
	if (!is.factor(class)) class<-as.factor(class)

	val.scale<-scale(val)	
	p<-dim(val)[2]			# Nb variables
	Maxi<-max(val.scale,na.rm=T)
	mini<-min(val.scale,na.rm=T)
	l<-length(levels(class))	# Nb niveaux dans la classe
	increm<-0.1*(p-1)			# incrémentation: décalage des boîtes


	plot(val.scale[,1]~class,boxwex=0.15,ylim=c(mini-1,Maxi+0.5),at=1:l-increm,xlab=xlab,ylab=ylab,
		names= substr(levels(class),start=1,stop=nchar),col=rainbow(100,1)[1],main=title)

	for (i in 2:p)
	{
		plot(val.scale[,i]~class,boxwex=0.15,ylim=c(mini-1,Maxi+0.5),at=1:l-increm+0.15*(i-1),add=TRUE,
		col=rainbow(10,1)[i],xlab=NULL,ylab=NULL,names=substr(levels(class),start=1,stop=nchar))
	}
	abline(h=0,col="blue",lty=3)

	if (p<=3)	legend("bottomleft",colnames(val.scale),fill=rainbow(10,1)[1:3],cex=0.75,bty="n")

	else	{	
		legend("bottomleft",colnames(val.scale)[1:3],fill=rainbow(10,1)[1:3],cex=0.75,bty="n")
		legend("bottomright",colnames(val.scale)[4:p],fill=rainbow(10,1)[4:p],cex=0.75,bty="n")
		}
}




