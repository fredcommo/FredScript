ROC<-function(class,val,nclass,title=NULL)

{
	if(length(which(is.na(val)))>0)	# suppression des  NA's
	{
	index.na<-which(is.na(val))	
	val<-val[-index.na]
	class<-class[-index.na]
	}
	
	
	seuil<-rep(0,(nclass))
	err<-rep(0,(nclass))		
	sens<-rep(0,(nclass))
	vpp<-rep(0,(nclass))
	vpn<-rep(0,(nclass))

	for (i in 1:(nclass))		# détermine les valeurs seuil
	seuil[i]<-max(val)/nclass*i
	

	for (i in 1:(nclass))		# réaffecte une valeur en fct du seuil
		{	
		tmp<-ifelse(val<seuil[i],"moins","plus")
		tabl<-table(tmp,class)
		err[i]<-tabl[2,1]/colSums(tabl)[1]		# calcule l'erreur: P(plus|nonM)
		sens[i]<-tabl[2,2]/colSums(tabl)[2]		# calcule le nbr de bien classés:P(plus|M)
		vpp[i]<-tabl[2,2]/rowSums(tabl)[2]		# valeur predictive positive:P(M|plus)
		vpn[i]<-tabl[1,1]/rowSums(tabl)[1]		# valeur predictive negative:P(nonM|moins)
		}

	resultat<-cbind(Seuil=seuil,Sensib=sens,Error=err)
	resultat<-as.data.frame(resultat)
	
	plot(resultat$Sensib~resultat$Error,type="l",
		xlim=range(0,1),ylim=range(0,1),
		xlab="1-spécificité",ylab="Sensibilité",
		col="blue",main=title)

	x<-seq(0,1,len=nclass)		# ajoute y=x
	y<-x
	lines(y~x)

		# Calcule les distances de la courbe à la droite y=x
	resultat<-cbind(resultat,x=seq(1,0,len=nclass),y=seq(1,0,len=nclass))
	distance<-sqrt((resultat$Sensib-resultat$x)^2+(resultat$Error-resultat$y)^2)
	resultat<-cbind(resultat,dist.XY=distance)

		# Trace les Sens et Erreur à la meilleure valeur
	max.dist<-max(resultat$dist.XY)
	best<-which(resultat$dist.XY==max.dist)

	if (length(best)>1)
		{index.best<-round(length(best)/2)
		best<-best[index.best]
		}
	
	best.sens<-resultat$Sensib[best]
	best.err<-resultat$Error[best]
	best.vpp<-vpp[best]
	best.vpn<-vpn[best]
	cut<-resultat$Seuil[best]

	# Légendes
	abline(h=best.sens,v=best.err,lty=3,col="red")
	legend(x=best.err,y=best.sens,legend=paste("value=",signif(cut,4),sep=""),cex=0.75,text.col="blue",bty="n")
	legend(x=-0.06,y=best.sens+0.05,legend=paste("Sens=",signif(best.sens,2)*100,"%",sep=""),cex=.7,text.col="red",bty="n")
	legend(x=best.err-0.03,y=0.03,legend=paste("Err=",signif(best.err,2)*100,"%",sep=""),cex=.7,text.col="red",bty="n")
	legend("bottomright",legend=c(paste("VPP=",signif(best.vpp,2),sep=""),paste("VPN=",signif(best.vpn,2),sep="")),cex=.7,text.col="red",bty="n")
	
	list(Tab=resultat[,-c(4:5)],cut=cut)					# renvoie la meilleure valeur
}