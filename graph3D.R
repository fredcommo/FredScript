graph3D<-function(x,y,z,xlab=NULL,ylab=NULL,zlab=NULL,expand.size=2)

# Représentation graphique en 3D, à partir de la fonction cloud(), du résultat d'une ACP

# ACP 		: Résultat d'une ACP (ACP<-prcomp(data))
# class1,class2 	: identifiant des groupes. De type factor ou numeric (les valeurs seront converties en facteurs dans ce cas). 
#			  Peur aussi être le résultat d'un kmeans() ou autre clustering.
# expand.size 	: facteur ajustable d'expansion de taille des points.

{
library(lattice)					# nécessaire pour cloud().


# Définit l'origine
	x0<-min(x)
	y0<-min(y)
	z0<-max(z)

# Calcul des distances à l'origine (x0,y0,z0)

	d1<-sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)
	d2<-log(max(d1)/d1)

	cloud(z~x*y,cex=expand.size*d2,
			xlab=xlab,ylab=ylab,zlab=zlab)

}