kinome.analysis.v2<-function(X,Y)

# X    	: reference experiments (untreated)
# Y	  	: treated datas

{
	ref<-cbind.data.frame(value1=X$value1,value2=X$value2)
	exp<-cbind.data.frame(value1=Y$value1,value2=Y$value2)
	wells<-X$well
	kinases<-X$kinase
	substrates<-X$substrate

	#p<-dim(ref)[2]
	
	lev<-length(wells)
	init<-rep(0,lev)	

	# Fichier de sortie
	result<-cbind.data.frame(well=wells,kinase=kinases,
						Substrate=substrates,
						mu.Ref=init, mu.Exp=init,Ratio=init,
						Log.R=init,foldchange=init,
						p.value=init,adj.p=init)
	
	# Affectation des "well" en tête de ligne
	#i=1
	#for (j in 1:16)
	#	for (k in 1:24)
	#	{
	#		result[i,1]<-paste(LETTERS[j],k,sep="")
	#		i=i+1
	#	}	

	# Remplacement des valeurs négatives (renommées NA) par 1 pour calcul Log
	ref[which(is.na(ref[,1])),1:2]<-1
	exp[which(is.na(exp[,1])),1:2]<-1

	# Option si neg values non supprimées
	ref[which(ref[,1]<=0),1]<-1
	ref[which(ref[,2]<=0),2]<-1
	exp[which(exp[,1]<=0),1]<-1
	exp[which(exp[,2]<=0),2]<-1

	Log.R<-rep(0,lev)

		for (i in 1:lev)
		{
			mu1<-mean(as.numeric(ref[i,]))
			mu2<-mean(as.numeric(exp[i,]))
			Log.R[i]<-round(log2(mu2/mu1),2)		# Calcul du Log.ratio (base=2)
	
			result$mu.Ref[i]<-mu1
			result$mu.Exp[i]<-mu2
			result$Ratio[i]<-signif(mu2/mu1,2)
		}
		
		mu.log<- mean(Log.R)
		sd.log<-sd(Log.R)

		# Calcul des p.values et p.values ajustées
		result$Log.R<-Log.R
		result$foldchange<-ifelse(result$Ratio<1,
							round(-1/result$Ratio,0),
							round(result$Ratio,0))

		result$p.value<-round(1-pnorm(abs(Log.R),mean=mu.log,sd=sd.log),4)	#*2
		result$adj.p<-round(p.adjust(result$p.value,method="BH"),4)
		
		result
}
