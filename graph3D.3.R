graph3D.3<-function(X,class1=NULL,class2=NULL,size=2, log.base=4, bold=TRUE,title=NULL)

#  3D-plot for ACP result.

# X 			: ACP result. X<-prcomp() or X<-princomp, but better results with prcomp()
# class1,class2 	: Class Ids. May be factor or numeric. In that case the values will be convert to factors. 
#			  May be the result of any clustering method : kmeans(),LDA(),cutree(hclust_product)...
# size 		: graphic parameter to expand the point.
# bold		: graphic parameter. Default is 'TRUE'.
# title		: graph title


# Begin
{
	# Package "lattice" contains cloud() for the 3D-plot.
		library(lattice)	
				
	# Convert classes to factors if they are not.
		if (is.numeric(class1))				
			class1<-as.factor(class1)
		if (is.numeric(class2))
			class2<-as.factor(class2)

		n1<-length(levels(class1))
		n2<-length(levels(class2))

	# Define axes and compute scaled values (proportion of max-values)

	if (class(X)=="prcomp")
		{
			x<-X$x[,1]/max(X$x[,1])
			y<-X$x[,2]/max(X$x[,2])
			z<-X$x[,3]/max(X$x[,3])
			name.x="PC1"; name.y="PC2"; name.z="PC3"
		}

	else
		{
			X<-as.matrix(X)
			x<-X[,1]/max(X[,1])
			y<-X[,2]/max(X[,2])
			z<-X[,3]/max(X[,3])
			name.x=colnames(X)[1]; name.y=colnames(X)[2]; name.z=colnames(X)[3]			
		}
	


	# Define the distance origin as the farest corner

		x0<-max(x)
		y0<-max(y)
		z0<-min(z)

	# Compute distances from the origin (x0,y0,z0)

		d1<-sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)


	# Define color values

		col.palette<-c("royalblue3","indianred2", "blue2", "firebrick1", "chartreuse3", "deeppink1","forestgreen","deepskyblue2","goldenrod2", "purple1", "burlywood3")	# 1 = "black" is excluded	
		if (!is.null(class1)) color=col.palette[1:n1][class1]
			
		else
			{
				n1=1
				color="blue2"
			}


	# Define points values

		if (bold) pts.palette<-c(16:20,15)		
		else pts.palette<-c(1,2,5,6,21,22)

		if (is.null(class2)) points=pts.palette[1:n1][class1]
		else points=pts.palette[1:n2][class2]
			

	# draw the 3D graph
		cloud(z~x*y,pch=points,col=color,cex=log(d1,base=log.base)*size,xlab=name.x,ylab=name.y,zlab=name.z, main=title)
		#cloud(z~x*y,pch=points,col=color,cex=d1*size,xlab=name.x,ylab=name.y,zlab=name.z, main=title)
				
			#xlim=range(0,max(x)),ylim=range(0,max(y)),zlim=range(0,max(z)),

} 

# END

# EXAMPLE 1: Iris data, 4 numeric variables (sizes) and 1 factor variable (3 different Species)

# 	source("D:\\Stats\\Doc R\\Scripts R\\graph3D.3.R") 			# source from script folder
# 	data(iris)
#	acp.iris<-prcomp(iris[,-5])							# column 5 contains species, so not numeric.
# 	graph3D.3(acp.iris, class1=iris$Species, bold=F)


# EXAMPLE 2: Crabs data, 5 numeric variables (sizes) and 2 factor variable (sex, sp)

#	library(MASS)	#(contains crabs datas)
#	graph3D.3(prcomp(crabs[,-c(1:3)]/crabs[,5]), class1=crabs$sp, class2=crabs$sex, size=5, bold=T)




