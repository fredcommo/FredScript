kinome.analysis<-function(X,Y)

# X    	: reference experiments (untreated)
# Y	  	: treated datas

{
	p<-dim(X)[2]
	lev<-levels(X$well)
	n.lev<-length(lev)-1

	# Fichier de sortie
	result<-cbind.data.frame(well=X$well[1],kinase=X$kinase[1],
						Substrate=X$substrate[1],
						mu.X=0, mu.Y=0,Ratio=0,
						Log.R=0,foldchange=0,p.value=0,adj.p=0)
	
	# Affectation des "well" en tête de ligne
	i=1
	for (j in 1:16)
		for (k in 1:24)
		{
			result[i,1]<-paste(LETTERS[j],k,sep="")
			i=i+1
		}	

	# Remplacement des valeurs négatives (en fait <=1) par 1 pour calcul Log
	X$value[which(X$value<=1)]<-1
	Y$value[which(Y$value<=1)]<-1
	
	Log.R<-rep(0,n.lev)

		for (i in 1:n.lev)
		{
			well<-result$well[i]
	
			mu1<-mean(X[which(X$well==well),]$value)
			mu1<-round(mu1,2)
			mu2<-mean(Y[which(Y$well==well),]$value)
			mu2<-round(mu2,2)
		
			Log.R[i]<-round(log2(mu2/mu1),2)		# Calcul du Log.ratio (base=2)
	
			result$kinase[i]<-X$kinase[which(X$well==well)][1]
			result$Substrate[i]<- X$substrate[which(X$well==well)][1]
			result$mu.X[i]<-mu1
			result$mu.Y[i]<-mu2
			result$Ratio[i]=signif(mu2/mu1,2)
		}
		
		mu.log<- mean(Log.R)
		sd.log<-sd(Log.R)

		# Calcul des p.values et p.values ajustées
		result$Log.R<-Log.R
		result$foldchange<-ifelse(result$Ratio<1,
							round(-1/result$Ratio,0),
							round(result$Ratio,0))

		result$p.value<-1-pnorm(abs(Log.R),mean=mu.log,sd=sd.log)
		result$adj.p<-p.adjust(result$p.value,method="BH")
		
		result
}
