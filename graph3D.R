graph3D<-function(x,y,z,xlab=NULL,ylab=NULL,zlab=NULL,expand.size=2)

# Repr�sentation graphique en 3D, � partir de la fonction cloud(), du r�sultat d'une ACP

# ACP 		: R�sultat d'une ACP (ACP<-prcomp(data))
# class1,class2 	: identifiant des groupes. De type factor ou numeric (les valeurs seront converties en facteurs dans ce cas). 
#			  Peur aussi �tre le r�sultat d'un kmeans() ou autre clustering.
# expand.size 	: facteur ajustable d'expansion de taille des points.

{
library(lattice)					# n�cessaire pour cloud().


# D�finit l'origine
	x0<-min(x)
	y0<-min(y)
	z0<-max(z)

# Calcul des distances � l'origine (x0,y0,z0)

	d1<-sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)
	d2<-log(max(d1)/d1)

	cloud(z~x*y,cex=expand.size*d2,
			xlab=xlab,ylab=ylab,zlab=zlab)

}