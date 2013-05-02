scurve2 <- function(x, y, n = 2, bw = 0.2, output = F, Col = "blue", LWD = 3){

f <- function(p, x){
	y <- p[1]*x^2 + p[2]*x + p[3]
	return(y)
	}

sce <- function(p, xobs, yobs){
	yfit <- f(p, xobs)
	res <- sum((yfit-yobs)^2)
	return(res)
	}

	init <- c(1, 1, 1)

	X <- c()
	Y <- c()

	for (i in 1:(length(x) - n)){
		first = i
		last = i + n

		opt <- optim(p = init, sce, xobs=x[first:last], yobs=y[first:last])

		newx <- seq(x[first], x[last], len=50)
		newy <- f(opt$par, newx)
		X <- c(X, newx)
		Y <- c(Y, newy)
		}		
	lines(ksmooth(X, Y, kernel = "normal", bandwidth = bw), col = Col, lwd = LWD)
	if (output) return(list(x = X, y = Y))
}
