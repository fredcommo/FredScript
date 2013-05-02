						##################

# A partir d'un fichier abr3.txt, la fonction renvoie un data.frame trié par band1
# La valeur Log2Ratio est la moyenne des Log2Ratios des régions (100Kb) couvrant cette bande.

						###################

# Possibilités:

	# hclust() par patient directement sur le fichier de sortie ou par région
	# hclust() par région sur le fichier de sortie transposé: t(fichier_sortie)
	# boxplot() par région sur le fichier de sortie

						####################


MeanbyBand<-function(x,...)				
{
	data<-x;

	fact<-unique(factor(data$band1));		# Renvoie les niveaux des facteurs Band1
	l<-length(unique(factor(data$band1)));

	names<-colnames(data[,7:23]);			# Renvoie les noms d'exp
	n<-length(names);

	tab<-matrix(0,n,l);				# Crée une matrice(n,l) pour recevoir les res.

	rownames(tab)<-names;				# Nomme les lignes et colonnes
	colnames(tab)<-fact;

	for (i in 1:l)					# Boucle de calcul des moyennes après extraction des bandes
		{
			
			tab[,i]<-mean(data[data$band1==fact[i],7:23]);
			
		}
	
	tab<-as.data.frame(tab);


}