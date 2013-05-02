Summary.Unique.Chr<-function(x,chrom.number)

# Calcule les fréquences de gain et perte de l'ensemble de l'effectif sur les sequences d'un chromosome donné.

	# x : objet de type aCGH
	# chrom.number : n° du chromosome

{
	aCGH.obj<-x

	index<-which(clones.info(aCGH.obj)$Chrom==chrom.number)	
	aCGH.Chr<-aCGH.obj[index, keep=T]

	summarize.clones(aCGH.Chr)->Chr.summary

	Chr.summary
}
