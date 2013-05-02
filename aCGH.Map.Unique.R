aCGH.Map.Unique<-function(x,chrom.selection)

# x : objet aCGH
# chrom.selection : #Chr à représenter

	{

		aCGH.Obj<-x
		Gain.Loss.Ch.XX<-data.frame(0)


	# Création d'un index (exemple: création d'une sous-liste relative au 6q) et extraction selon cet index.
	# Le paramètre "keep=T" permet d'associer les éléments de la liste principale dans la sous-liste aCGH.Ch.XX
	
		index<-which(clones.info(aCGH.Obj)$Chrom==chrom.selection)	
		aCGH.Ch.XX<-aCGH.Obj[index, keep=T]
	
		aCGH.Ch.XX$clones.info$Chrom<-rep(1,length(clones.info(aCGH.Ch.XX)$Chrom))

		plotFreqStat(aCGH.Ch.XX, X=F, Y=F, titles=paste("Chrom ",chrom.selection), numaut=1, chrominfo=human.chrom.info.Jul03[chrom.selection,])

	#abline(h=freq1,lty=3,col="blue")	# ajout de lignes horizontales supplémentaires (repères de fréquences)
	#abline(h=freq2,lty=3,col="red")


	}			