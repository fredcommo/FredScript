plot.replic<-function(variable,resp,xlab=NULL,title=NULL)

# Représente les valeurs moyennes de réponse et sdev par niveau de variable (réplicats)

{
	lev<-unique(variable)
	nval<-length(lev)

	mean.val<-rep(0,nval)
	sdev<-rep(0,nval)

	# Calcul des moyennes et sd pour chaque niveau
	for (i in 1:nval)
		{	mean.val[i]<-mean(resp[variable==lev[i]])
			sdev[i]<-sd(resp[variable==lev[i]])
		}
	
	# Représente les valeurs moyennes par niveau
	plot(lev,mean.val,xlab=xlab,ylab="mean resp",main=title)
	
	# Ajoute les barres d'erreur (+/-sd)

	size<-abs(max(lev)-min(lev))/100					# mise à l'échelle des barres horizontales

	segments(lev,mean.val-sdev,lev,mean.val+sdev)			# barres verticales
	segments(lev-size,mean.val-sdev,lev+size,mean.val-sdev)	# barres horizontales basses
	segments(lev-size,mean.val+sdev,lev+size,mean.val+sdev)	# barres horizontales hautes

	# Renvoie les valeurs moyennes et sd
	cbind.data.frame(values=lev,mean=mean.val,sdev=sdev)
}