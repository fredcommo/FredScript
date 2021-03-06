MultiBox2.scale<-function(class, val, xlab="class", ylab="values", nchar=10, title=NULL, legcex=0.75)

# fonction graphique pour la repr�sentation boxplot de plusieurs variables
# en fonction des niveaux de classe

# A partir des donn�es centr�es r�duites

{
	if (!is.factor(class)) class<-as.factor(class)
	
	val.scale<-scale(val)	
	p<-dim(val)[2]			# Nb variables
	p1<-round(p/2,0)
	x<-p*length(levels(class))
	wex<-1/x
	Maxi<-max(val.scale,na.rm=T)
	mini<-min(val.scale,na.rm=T)
	l<-length(levels(class))	# Nb niveaux dans la classe
	increm<-0.01*(p-1)			# incr�mentation: d�calage des bo�tes


	plot(val.scale[,1]~class,boxwex=wex,ylim=c(mini-1,Maxi+0.5),at=1:l-increm,axes=FALSE, xlab=xlab,
ylab=ylab,names= substr(levels(class),start=1,stop=nchar),col=rainbow(p,1)[1],main=title)

	for (i in 2:p1)
	{
plot(val.scale[,i]~class,boxwex=wex,ylim=c(mini-1,Maxi+0.5),at=1:l-increm+wex*(i-1),add=TRUE, axes=FALSE, 
col=rainbow(p,1)[i],xlab=NULL,ylab=NULL,names=substr(levels(class),start=1,stop=nchar))
	}

plot(val.scale[,(p1+1)]~class,boxwex=wex,ylim=c(mini-1,Maxi+0.5),
at=1:l-increm+wex*p1,add=TRUE, col=rainbow(p,1)[(p1+1)],xlab=NULL,ylab=NULL,
names=substr(levels(class),start=1,stop=nchar))

	for (i in (p1+2):p)
	{
plot(val.scale[,i]~class,boxwex=wex,ylim=c(mini-1,Maxi+0.5),at=1:l-increm+wex*(i-1),add=TRUE, axes=FALSE,
col=rainbow(p,1)[i],xlab=NULL,ylab=NULL,names=substr(levels(class),start=1,stop=nchar))
	}


	abline(h=0,col="blue",lty=3)

	if (p<=3)	
legend("bottomleft",colnames(val.scale),fill=rainbow(p,1)[1:3],cex=legcex,bty="n")

	else	{	
legend("bottomleft",colnames(val.scale)[1:p1],fill=rainbow(p,1)[1:p1],cex=legcex, bty="n")
legend("bottomright",colnames(val.scale)[c(p1+1):p],fill=rainbow(p,1)[(p1+1):p],cex=legcex, bty="n")
		}

}




