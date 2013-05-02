BoxPoints <- function(x, y, densit=FALSE, Med = F, disp=0.05,...){

# x : Class (as factor)
# y : values
# densit : if true, a density line is added on the plot (vertically)
# disp : a graphic factor to set the dispersion of points.

require(affy)
# library(affy) for tukey.biweight()

	x <- as.factor(x)
	if(any(c(is.na(x), is.na(y)))){
		na.index <- which(is.na(x) | is.na(y))
	 	x <- x[-na.index]
		y <- y[-na.index]
		}

	tmp <- cbind.data.frame(x,y)
	colnames(tmp) <- c("x","y")
	tmp <- tmp[order(x),]
	x <- tmp$x
	y <- tmp$y

	X <- c()
	med <- rep(0, nlevels(x))

	for (i in 1:nlevels(x)){
		X <- c(X, rnorm(length(x[x==levels(x)[i]]),i, disp))
		med[i] <- tukey.biweight(y[x==levels(x)[i]])
	}

	leg.pos <- "topleft"
	nlev <- nlevels(x)
	if(max(y[x==levels(x)[1]])>max(y[x==levels(x)[nlev]])) leg.pos <- "topright"

	col.palette = rainbow(length(y), start = 0.6, end = 0.10)
	pt.col= as.numeric(rank(y))

	# if(is.na(Ylim)) Ylim=range(min(y),max(y))

	plot(y ~ x, border = NA,...)
	points(y ~ X, col=col.palette[pt.col],...)	# pch=c(2:(nlevels(x)+1))[x]

	if(Med){
		segments(x0=seq(1:nlevels(x))-0.15, y0=med, x1=seq(1:nlevels(x))+0.15, y1=med, lwd=3, col="violetred4")
		legend(leg.pos, legend = "medians", lwd = 3, col = "violetred4", bty = "n", cex=1)
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