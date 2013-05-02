graph3D.9<-function(X, class1 = NULL, class2 = NULL, size = 0.5, Cols = NULL, log.base = 2, bold = TRUE, title = NULL, add.legend = T, Theta = 0, Phi = 0, Z = 0,...)

#  3D-plot for ACP result.

# X 			: ACP result. X<-prcomp() or X <- princomp, but better results with prcomp()
# class1,class2 	: Class Ids. May be factor or numeric. In that case the values will be convert to factors. 
#			  May be the result of any clustering method : kmeans(),LDA(),cutree(hclust_product)...
# size 		: graphic parameter to expand the point.
# Cols		: define the colors for points according to factors. If 'NULL' (default) colors will be automatically defined by 'class1' levels, if !null.
# bold		: graphic parameter. Default is 'TRUE'.
# title		: graph title
# Theta, Phi, Z	: rotation angles (in radius), resp on X, Y and Z axis


# Begin
{
	# Package "lattice" contains cloud() for the 3D-plot.
		library(lattice)	

	# Build the rotation matrix
	Rx = matrix(c(1, 0, 0, 0, cos(Theta), sin(Theta), 0, -sin(Theta), cos(Theta)), 3, 3)
	Ry = matrix(c(cos(Phi), 0, -sin(Phi), 0, 1, 0, sin(Phi), 0, cos(Phi)), 3, 3)	
	Rz = matrix(c(cos(Z), -sin(Z), 0, sin(Z), cos(Z), 0, 0, 0, 1), 3, 3)	
	R = Rx%*%Ry%*%Rz

	# Define axes and compute the rotation values
	{
	if (class(X)=="prcomp") {
		if(ncol(X$x)<3) stop(paste("Your matrix has only", ncol(X$x), "columns. At least 3 are expected !\n"))
		if(ncol(X$x)>3) cat(paste("Your matrix has", ncol(X$x), "columns. The 3 firsts have been considered !\n"))
		name.x = colnames(X$x)[1]; name.y = colnames(X$x)[2]; name.z = colnames(X$x)[3]
		X <- as.matrix(X$x[,1:3])%*%R
		}
	else {
		if(ncol(X)<3) stop(paste("Your matrix has only", ncol(X), "columns. At least 3 are expected !\n"))
		if(ncol(X)>3) cat(paste("Your matrix has", ncol(X$x), "columns. The 3 firsts have been considered !\n"))
		name.x = colnames(X)[1]; name.y = colnames(X)[2]; name.z = colnames(X)[3]
		X <- as.matrix(X[,1:3])%*%R
		}
	}	
	x <- X[,1]
	y <- X[,2]
	z <- X[,3]

	# Convert classes to factors if they are not.
		if (is.null(class1)) class1 <- rep(1, length(x))				
		if (is.null(class2)) class2 <- rep(1, length(x))

		class1 <- as.factor(class1)
		class2 <- as.factor(class2)

		n1 <- nlevels(class1)
		n2 <- nlevels(class2)



	# Define the distance origin as the farest corner

		x0 <- max(x, na.rm = T)
		y0 <- max(y, na.rm = T)
		z0 <- min(z, na.rm = T)

	# Compute distances from the origin (x0,y0,z0)

		d1 <- sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)


	# Define color values
		if(is.null(Cols)){
		Cols <- c("royalblue2", "indianred1", "seagreen3", "goldenrod2", "blue3", "purple1", "burlywood2", "firebrick1", "green4", "deeppink1", "deepskyblue2")	# 1 = "black" is excluded	
		}
		color = Cols[1:n1][class1]			


	# Define points values
		{
		if (bold) pts.palette <- c(16:20,15)		
		else pts.palette <- c(1,2,5,6,21,22)
		}	
		points = pts.palette[1:n2][class2]
	

	# draw the 3D graph
		x <- x[order(d1, decreasing = F)]
		y <- y[order(d1, decreasing = F)]
		color <- color[order(d1, decreasing = F)]
		points <- points[order(d1, decreasing = F)]
		z <- z[order(d1, decreasing = F)]
		d1 <- d1[order(d1, decreasing = F)]

		if(Theta+Phi+Z == 0) cloud(z~x*y, pch = points, col = color, cex = -log(d1/(max(d1)*1.05), base = log.base)*size, 
						xlab = name.x, ylab = name.y, zlab = name.z, main = title,...)
	
		else 			cloud(z~x*y, pch = points, col = color, cex = log(d1, base = log.base)*size, 
						xlab = NULL, ylab = NULL, zlab = NULL, default.scales = list(draw = F), main = title,...)


		# Legend <- NA
		# if(add.legend & (nlevels(class1)!=1 | nlevels(class2)!=1)){
		#		{
		#		if(levels(class2)==1) Legend <- levels(class1)
		#		else Legend <- levels(class1:class2)
		#		}
		#	if(any(!is.na(levels(Legend)))){
		#		points(x = seq(0.075, 0.975, len = nlevels(Legend)), y = rep(0, nlevels(Legend)), pch = points)
		#		text(x = seq(0.1, 1, len = nlevels(Legend)), y = rep(0, nlevels(Legend)), labels = levels(Legend))
		#		}
		#	}		
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




