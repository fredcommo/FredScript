						##################

# A partir d'un fichier abr3.txt, la fonction renvoie un data.frame tri� par band1
# La valeur Log2Ratio est la moyenne des Log2Ratios des r�gions (100Kb) couvrant cette bande.

						###################

# Possibilit�s:

	# hclust() par patient directement sur le fichier de sortie ou par r�gion
	# hclust() par r�gion sur le fichier de sortie transpos�: t(fichier_sortie)
	# boxplot() par r�gion sur le fichier de sortie

						####################


MeanbyBand<-function(x,...)				
{
	data<-x;

	fact<-unique(factor(data$band1));		# Renvoie les niveaux des facteurs Band1
	l<-length(unique(factor(data$band1)));

	names<-colnames(data[,7:23]);			# Renvoie les noms d'exp
	n<-length(names);

	tab<-matrix(0,n,l);				# Cr�e une matrice(n,l) pour recevoir les res.

	rownames(tab)<-names;				# Nomme les lignes et colonnes
	colnames(tab)<-fact;

	for (i in 1:l)					# Boucle de calcul des moyennes apr�s extraction des bandes
		{
			
			tab[,i]<-mean(data[data$band1==fact[i],7:23]);
			
		}
	
	tab<-as.data.frame(tab);


}