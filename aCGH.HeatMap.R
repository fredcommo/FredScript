aCGH.HeatMap<-function(x,cut=1,title="")

# Représentation en heatmap

	# x : objet de classe aCGH
	# title : titre du graphique

	{
	aCGH.Obj<-x
	clusterGenome(aCGH.Obj, dendPlot = F, titles = title, methodS = "ward", ncolors = 100,cutoff = cut, imp = F, 
	vecchrom = 1:22,samplenames = phenotype(aCGH.Obj)$Id)
	}
