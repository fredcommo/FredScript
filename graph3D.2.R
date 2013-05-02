graph3D.2<-function(X,class1=NULL,class2=NULL,size=1.5,bold=TRUE,title=NULL)

# Représentation graphique en 3D, à partir de la fonction cloud(), du résultat d'une ACP

# X 			: Résultat d'une ACP, X<-prcomp(data)
# class1,class2 	: identifiant des groupes. De type factor ou numeric (les valeurs seront converties en facteurs dans ce cas). 
#			  Peur aussi être le résultat d'un kmeans() ou autre clustering.
# size 		: paramètre graphique : facteur ajustable d'expansion de taille des points.
# bold		: paramètre graphique : si "True" points pleins, sinon vides

# Begin
{
	# nécessaire pour cloud().
		library(lattice)	
				
	# conversion des identifiants de groupes en facteurs si nécessaire.
		if (is.numeric(class1))				
			class1<-as.factor(class1)
		if (is.numeric(class2))
			class2<-as.factor(class2)

		n1<-length(levels(class1))
		n2<-length(levels(class2))

	# Définit les axes avec echelle proportionnelle
	{
	if (class(X)=="prcomp")
		{
		x<-X$x[,1]/max(X$x[,1])
		y<-X$x[,2]/max(X$x[,2])
		z<-X$x[,3]/max(X$x[,3])
		name.x="PC1"; name.y="PC2"; name.z="PC3"
		}
	else
		{
		x<-X[,1]/max(X[,1])
		y<-X[,2]/max(X[,2])
		z<-X[,3]/max(X[,3])
		name.x=colnames(X)[1]; name.y=colnames(X)[2]; name.z=colnames(X)[3]
		}
	}
	# Définit l'origine
		x0<-min(x)
		y0<-min(y)
		z0<-max(z)

	# Calcul des distances à l'origine (x0,y0,z0)

		d1<-sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)
		d2<-(-log(d1*0.10,base=10))

	# Définition des paramètres graphiques : couleurs
		col.palette<-c(1:n1)
			{
			if (is.null(class2) & !is.null(class1))
				color=col.palette[1:n1][class1]
			else  color=col.palette[1:n2][class2]
			}
			if (is.null(class1))
			{
				n1=1
				color="blue"
			}

	# Définition des paramètres graphiques : points
			if (bold) pts.palette<-c(16:20,15)		
			else pts.palette<-c(1,2,5,6,21,22)
		
			if (is.null(class2)) points=ifelse(is.null(class1),8,pts.palette[1:n1][class1])
			else points=pts.palette[1:n2][class2]
			

	# Trace le graphe
		cloud(z~x*y,pch=points,col=color,cex=size*d2,
				xlim=range(0,max(x)),ylim=range(0,max(y)),zlim=range(0,max(z)),
				xlab=name.x,ylab=name.y,zlab=name.z, main=title)

} 
# END

# data(iris)
# source("D:\\Stats\\Doc R\\Scripts R\\graph3D.2.R")
# prcomp(iris[,-5])->acp.iris
# graph3D.2(acp.iris, class1=iris$Species)


