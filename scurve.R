scurve <- function(x, y, n=3){

	for (i in 1:n){
	m <- length(x)-1
	newx <- c()
	newy <- c()
		for (j in 1:m){
		newx <- c(newx, x[j], (x[j]+x[j+1])/2)
		newy <- c(newy, y[j], (y[j]+y[j+1])/2)
		}
	x <- c(newx, x[length(x)])
	y <- c(newy, y[length(y)])
	}
	return(list(x=x, y=y))
}
