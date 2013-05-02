cellNames <- info$Characteristics.Cell.Line.Name.

k=1
for(i in seq(0,2,len=16)){
	plot(acp$x, type="n", xlab="", ylab="", xlim=range(-2.5, 3), ylim=range(-2, 2))
	t = i*pi

	Xp <- cbind(pc1, pc2, pc3)
	Xt <- cbind(rot1, rot2, rot3)
	Rx = matrix(c(1, 0, 0, 0, cos(t), sin(t), 0, -sin(t), cos(t)), 3, 3)		
	Ry = matrix(c(cos(t), 0, -sin(t), 0, 1, 0, sin(t), 0, cos(t)), 3, 3)		
	Rz = matrix(c(cos(t), sin(t), 0, -sin(t), cos(t), 0, 0, 0, 1), 3, 3)		
	Xp_rot = Xp%*%Ry
	Xt_rot = Xt%*%Ry
	mp = min(Xp_rot[,3])
	Mp = max(Xp_rot[,3])
	par(bg="transparent")
	# points(Xp_rot[,1], Xp_rot[,2], pch=19, col="royalblue", cex = 1.5*(Xp_rot[,3]-mp)/(Mp - mp) + 0.5)
	text(Xp_rot[,1], Xp_rot[,2], labels=cellNames, col="blue", cex = (Xp_rot[,3]-mp)/(Mp - mp) + 0.25)
	
	mt = min(Xt_rot[,3])
	Mt = max(Xt_rot[,3])
	len = 0.1*(Xt_rot[,3]-mt)/(Mt - mt)
	ang = 25*(Xt_rot[,3]-mt)/(Mt - mt)
	size = (Xt_rot[,3]-mt)/(Mt - mt) + 0.5
	coef = 2

	arrows(x0=0, y0=0, x1 = Xt_rot[,1]*coef, y1 = Xt_rot[,2]*coef, length = len, angle = ang, col="firebrick")
	text(x=Xt_rot[,1]*coef, y=Xt_rot[,2]*coef, labels=gettext(Symb[index]), font=4, cex= size, col="firebrick")

	savePlot(filename=paste("Image_", k, sep=""), type="tiff")
	k=k+1

}

