BoxPoints2 <- function(X, disp = 0.05,...){

# X : a (n,m) matrix (or data frame). The m columns will be vizualized.
# disp : a graphic factor to set the dispersion of points.


	n <- nrow(X)
	m <- ncol(X)

	# X <- log10(values)
	y <- c()
	for(i in 1:m) y <- c(y, X[,i])
	rk <- matrix(rank(y, ties.method = "min"), n, m)
	ncol <- length(unique(as.vector(rk)))
	col.palette <- rainbow(ncol, start = 0.7, end = 0.13)[factor(as.vector(rk))]
	col.palette <- matrix(col.palette, n, m)

	boxplot(X, col = "transparent", border = "transparent",...)
	for(i in 1:m) points(X[,i]~rnorm(n, i, disp), col = col.palette[,i],...)
	}
