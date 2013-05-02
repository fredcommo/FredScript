aCGH.Dend<-function(x,title="")

# Représentation en dendrogramme

	# x : objet de classe aCGH
	# title : titre du graphique

	{
	library(stats)
	
	aCGH.Obj<-x

	# Crée un vecteur Id, transpose la matrice log2.ratios et renomme les lignes

		as.vector(phenotype(aCGH.Obj)$Id)->Exp
		log2.ratios(x)->log2
		as.data.frame(t(log2))->transp
		Exp->rownames(transp)


		hclust(dist(transp),method="ward")->aCGH.cluster
		plot(aCGH.cluster,main=title)

	aCGH.cluster

	}
