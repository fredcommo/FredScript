ROC2<-function(class,val,nclass=500,D=NA,n=NA,title=NULL)

# class 	: defines the classes which have to be compared
# val 	: the values to test
# nclass	: a value to define the interspaces (thesholds) between tested values

# D		: which class represent the desease status	A rajouter dans le code. Attention à l'ordre des classes pour le calcul des probas
# n		: which class represent the normal status		A rajouter dans le code. Attention à l'ordre des classes pour le calcul des probas
# title	: the title for the graph (optional)

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

	Max.val<-max(val)
	min.val<-min(val)
	pas<-(Max.val-min.val)/nclass

	for (i in 1:(nclass))		# détermine les valeurs seuil
	seuil[i]<-min.val+pas*i
	
	# Attention à l'ordre des colonnes 
	for (i in 1:(nclass))		# réaffecte une valeur en fct du seuil
		{	
		tmp<-ifelse(val<seuil[i],"moins","plus")
		tabl<-table(tmp,class)
		sens[i]<-tabl[2,1]/colSums(tabl)[1]		# calcule le nbr de bien classés:P(plus|M)
		err[i]<-tabl[2,2]/colSums(tabl)[2]		# calcule l'erreur: P(plus|nonM)
		vpp[i]<-tabl[2,1]/rowSums(tabl)[2]		# valeur predictive positive:P(M|plus)
		vpn[i]<-tabl[1,2]/rowSums(tabl)[1]		# valeur predictive negative:P(nonM|moins)
		}

	resultat<-cbind(Seuil=seuil,Sensib=sens,Error=err)
	resultat<-as.data.frame(resultat)



	x<-seq(0,1,len=nclass)		# ajoute y=x
	y<-x
	plot(y~x, type="l", xlab="1-spécificité", ylab="Sensibilité")
	lines(resultat$Sensib~resultat$Error, col="blue", main=title)


		# Calcule les distances de la courbe au point (0,1)
	distance<-sqrt((1-resultat$Sensib)^2+(resultat$Error)^2)	#sqrt((new.x)^2+(1-new.y)^2)
	resultat<-cbind(resultat,dist=distance)

		# Trace les Sens et Erreur à la meilleure valeur
	min.dist<-min(distance)
	best<-which(distance==min.dist)
	
	if (length(best)>1)
		{index.best<-round(length(best)/2)
		best<-best[index.best]
		}
	
	best.sens<-resultat$Sensib[best]
	best.err<-resultat$Error[best]
	best.vpp<-vpp[best]
	best.vpn<-vpn[best]
	cut<-resultat$Seuil[best]
	
		# Légend best cut-off
	abline(h=best.sens,v=best.err,lty=3,col="red")
	legend(x=best.err,y=best.sens,legend=paste("value=",signif(cut,4),sep=""),cex=0.75,text.col="blue",bty="n")
	legend(x=-0.06,y=best.sens+0.05,legend=paste("Sens=",signif(best.sens,2)*100,"%",sep=""),cex=.7,text.col="red",bty="n")
	
		# VPP et VPN
	legend(x=best.err-0.03,y=0.03,legend=paste("Err=",signif(best.err,2)*100,"%",sep=""),cex=.75,text.col="red",bty="n")
	legend("bottomright",legend=c(paste("VPP=",signif(best.vpp,2),sep=""),paste("VPN=",signif(best.vpn,2),sep="")),cex=.75,text.col="red",bty="n")
	
	# table des résultats
	Tab<-table(val>cut,class)
	rownames(Tab)<-c(paste("inf",round(cut,0)),paste("Sup",round(cut,0)))

	list(Results=resultat,cut=cut,Table=Tab)			# renvoie la meilleure valeur
}
