flames <- function(xc, yc, n=100, dispx=0.1, dispy=0.5){
	y <- rnorm(n, yc, dispy)
	x <- rnorm(n, xc, dispx)

	col.palette = rainbow(length(y), start = 0.7, end = 0.15)
	pt.col= as.numeric(rank(y))

	points(x, y, pch=19, col=col.palette[pt.col])
}