aCGH.fga.Unique<-function(x,chrom.selection)

	# fichier de classe aCGH
	# n° du chromosome


{
	aCGH.Obj<-x
	Gain.Loss.Ch.XX<-data.frame(0)


# Création d'un index (exemple: création d'une sous-liste relative au 6q) et extraction selon cet index.
# Le paramètre "keep=T" permet d'associer les éléments de la liste principale dans la sous-liste aCGH.Ch.XX

	index<-which(clones.info(aCGH.Obj)$Chrom==chrom.selection)	
	aCGH.Ch.XX<-aCGH.Obj[index, keep=T]


	fga.func(aCGH.Ch.XX)->Ch.XX.fga
	PGain<-round(Ch.XX.fga$gainP*100,2)
	PLoss<-round(Ch.XX.fga$lossP*100,2)
	PTotal<-PGain+PLoss
	cbind(as.vector(phenotype(aCGH.Ch.XX)$Id),Gain.Prop=PGain,Loss.Prop=PLoss,Total.Prop=PTotal)->Gain.Loss.Ch.XX
	colnames(Gain.Loss.Ch.XX)[1]<-"Samples"

	as.data.frame(Gain.Loss.Ch.XX)

}