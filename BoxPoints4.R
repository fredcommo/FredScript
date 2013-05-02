BoxPoints4 <- function(x, y, densit=FALSE, add.stat = c("none", "add.med", "add.mean"), add.sem = T, ColType = c("C", "G"), disp=0.05, pt.size=1, title=NULL, Xlab = "", Ylab = "Values", Ylim=NA){

# x : Class (as factor)
# y : values
# densit : if true, a density line is added on the plot (vertically)
# disp : a graphic factor to set the dispersion of points.

require(affy) 	# for tukey.biweight()

	x <- as.factor(x)
	if(any(c(is.na(x), is.na(y)))){
		na.index <- which(is.na(x) | is.na(y))
	 	x <- x[-na.index]
		y <- y[-na.index]
		}
	n <- length(y)
	m <- nlevels(x)

	tmp <- cbind.data.frame(x,y)
	colnames(tmp) <- c("x","y")
	tmp <- tmp[order(x),]
	x <- tmp$x
	y <- tmp$y
	rk <- rank(y, ties.method = "min")


	X <- c()
	med <- rep(0, nlevels(x))

	for (i in 1:nlevels(x)){
		X <- c(X, rnorm(length(x[x==levels(x)[i]]),i, disp))
		med[i] <- tukey.biweight(y[x==levels(x)[i]])
	}

	leg.pos <- "topleft"
	if(max(y[x==levels(x)[1]])>max(y[x==levels(x)[m]])) leg.pos <- "topright"

	col.type <- match.arg(ColType)
	switch(col.type,  C = (col.palette = rainbow(n, start = 0.7, end = 0.13)[factor(rk)]),
				G = (col.palette = grey(seq(1/m, 1-1/(2*m), len = m))[x]))
	# switch(col.type,  C = (ByRow = F), G = (ByRow = T))
	
	# col.palette <- matrix(col.palette, n, m, byrow = ByRow)

	if(is.na(Ylim)) Ylim = range(min(y),max(y))

	plot(y ~ x, border = NA, xlab = Xlab, ylab = Ylab, ylim = Ylim, main = title)
	points(y ~ X, col = col.palette, pch=19, cex=pt.size)	# pch=c(2:(nlevels(x)+1))[x]

	stat.col = "black"
	Stat <- match.arg(add.stat)
	if(Stat!="none"){
	switch(Stat, add.med = (stat.fun = tukey.biweight), add.mean = (stat.fun = mean))
	switch(Stat, add.med = (stat.name = "medians"), add.mean = (stat.name = "means"))
		for (i in 1:m){
			stat <- stat.fun(y[x==levels(x)[i]])
			segments(x0 = i-0.1, y0 = stat, x1 = i+0.1, y1 = stat, lwd = 2, col = stat.col)
			}
	legend(leg.pos, legend = stat.name, lwd = 3, col = stat.col, bty = "n", cex = 1)
	}

	if(Stat=="none") add.sem = F
	if(add.sem){
	sem.col = "black"; lwd.seg = 2
	for(i in 1:m){
		xbar <- mean(y[x == levels(x)[i]])
		Sdev <- sd(y[x == levels(x)[i]])
		n <- length(which(x == levels(x)[i]))
		sem <- Sdev/sqrt(n)
		t <- qt(0.975, n-1)
		segments(x0 = i-0.05, y0 = xbar-t*sem, x1 = i+0.05, y1 = xbar-t*sem, lwd = lwd.seg, col = sem.col)
		segments(x0 = i-0.05, y0 = xbar+t*sem, x1 = i+0.05, y1 = xbar+t*sem, lwd = lwd.seg, col = sem.col)
		segments(x0 = i, y0 = xbar-t*sem, x1 = i, y1 = xbar+t*sem, lwd = lwd.seg, col = sem.col)
		}
	}

	if (densit){
		fact=nlevels(x)
		if (nlevels(x)==1) fact=1.5
		d<-density(y, na.rm=TRUE, n=length(y))
		dx<- d$x
		dy <- d$y/max(d$y)*fact
		lines(dx~dy, col="blue")	
	}
}	