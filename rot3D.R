# 3D rotation

# X : A matrix(x, y, z)
# rotation : choice of rotation axis
# Col : Vector of colors; Same length as nrow(X)
# ptSize : a factor to expand the size of points
# n.image : Number of images to generate.
# folder : the filepath

#	Rx = matrix(c(1, 0, 0, 0, cos(t), sin(t), 0, -sin(t), cos(t)), 3, 3)		
#	Ry = matrix(c(cos(t), 0, -sin(t), 0, 1, 0, sin(t), 0, cos(t)), 3, 3)		
#	Rz = matrix(c(cos(t), sin(t), 0, -sin(t), cos(t), 0, 0, 0, 1), 3, 3)		



rot3D <- function(X, rotation = c("x.axis", "y.axis", "z.axis"), Col = rep("blue", nrow(X)), ptSize = 1.5, Asp = NA,  n.image=16, folder = getwd()){


	k=1
	for(i in seq(0, 2, len=n.image)){

		plot(X[,1:2], type="n", bty="n", axes = F, xlab = "", ylab = "", asp = Asp) #, xlim=range(-2.5, 3), ylim=range(-2, 2))
		t = i*pi

		Xp <- as.matrix(X[,1:3])
		
		type <- match.arg(rotation)
		switch(type,
		x.axis = (R = matrix(c(1, 0, 0, 0, cos(t), sin(t), 0, -sin(t), cos(t)), 3, 3)),		
		y.axis = (R = matrix(c(cos(t), 0, -sin(t), 0, 1, 0, sin(t), 0, cos(t)), 3, 3)),		
		z.axis = (R = matrix(c(cos(t), sin(t), 0, -sin(t), cos(t), 0, 0, 0, 1), 3, 3)))		
			
		Xp_rot = Xp%*%R
		mp = min(Xp_rot[,3])
		Mp = max(Xp_rot[,3])

		pCol <- Col[order(Xp_rot[,3], decreasing = F)]
		Xp_rot <- Xp_rot[order(Xp_rot[,3], decreasing = F),]

		par(bg="transparent")
		# for (p in 1:nrow(Xp_rot))
			points(Xp_rot[,1], Xp_rot[,2], pch = 19, col = pCol, cex = ptSize*(Xp_rot[,3] - mp)/(Mp - mp) + 0.5)

		# text(Xp_rot[,1], Xp_rot[,2], labels=cellNames, col="blue", cex = (Xp_rot[,3]-mp)/(Mp - mp) + 0.25)
	
		savePlot(filename=paste(folder, "/Image_", k, sep=""), type="tiff")
		k=k+1
	}	
		par(bg = "white")
}
