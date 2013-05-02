plot.levels<- function(x, y, col=col, plot.type="Plot", xlab=NA, ylab=NA){

	x <- as.numeric(x)
	y <- as.numeric(y)

	p <- length(unique(x))

	m <- s <- rep(0, p)
	
	for (i in 1:p){
		m[i] <- mean(y[x==unique(x)[i]], na.rm=T)
		s[i] <- sd(y[x==unique(x)[i]], na.rm=T)
	}

	if (plot.type=="Plot") plot(m~unique(x), type="l", col=col, xlab=xlab, ylab=ylab)
	else lines(m~unique(x), col=col, xlab=xlab, ylab=ylab)

	segments(unique(x)-0.01, m-s, unique(x)+0.01, m-s)
	segments(unique(x)-0.01, m+s, unique(x)+0.01, m+s)
	segments(unique(x), m-s, unique(x), m+s)
}