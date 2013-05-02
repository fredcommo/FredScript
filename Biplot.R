Biplot<-function(x, choice = c(1:2), ptag = NULL, pt.cex = 2, pfont = 4, arrowlab= NULL, arrowcol = "darkred", exp.arrows = 1.5, Len = 0.15, Ang = 25,...){

	xdata <- x$x[,choice]
	xrot <- x$rotation[,choice]
	
	pcex = pt.cex
	if(ncol(xdata)>2){
		PC3 <- x$x[,3]
		mp = min(PC3)
		Mp = max(PC3)
		pcex = pt.cex*(PC3-mp)/(Mp - mp) + 0.5
		}

	mrot <- max(xrot)-min(xrot)
	mdata <- max(xdata[,choice])-min(xdata[,choice])

	coef <- mdata/mrot
	xrot <- xrot*coef/2

	choice1 <- choice[1]
	choice2 <- choice[2]

	x.range <- range(xdata[,choice1], xrot[,choice1])*1.05
	y.range <- range(xdata[,choice2], xrot[,choice2])*1.05
	
	palette <- c("goldenrod3", "navy", "indianred3", "seagreen3", "purple2", "paleturquoise4", "royalblue3")
	if (is.null(ptag)) {plab = 19; pcol = "royalblue3"}
	else {ptag <- as.factor(ptag); plab = ""; pcol = palette[ptag]}

	plot(xdata[,choice], pch = plab, col= pcol, cex = pcex, xlim = x.range, ylim = y.range,...)
	if(!is.null(ptag))
		text(x = xdata[,choice1], y = xdata[,choice2], labels = as.character(ptag) , col = pcol, cex = pcex, font = pfont)

	arrows(x0 = rep(0, nrow(xrot)), y0 = rep(0, nrow(xrot)),
		x1 = xrot[,choice1]*exp.arrows, y1 = xrot[,choice2]*exp.arrows, col = arrowcol, length = Len, angle = Ang,...)
	text(x = xrot[,choice1]*1.05, y = xrot[,choice2]*1.05, labels = arrowlab , col = arrowcol, font = pfont)
}