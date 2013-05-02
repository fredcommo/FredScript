carteKoho<-function(x, n, m,...)

{
	require(class)

	data <- x
	sg <- somgrid(n, m)
	rad <- c(seq(4, 0, len=10), rep(0,5))

	bSom <- batchSOM(data, sg, rad)

	class <- knn1(bSom$codes, data, 1:(n*m))
	bins <- as.numeric(class)

	plot(bSom$grid, type="n", xlim = range(0, m+1))
	symbols(bSom$grid$pts[,1], bSom$grid$pts[,2], circles=rep(0.5,n*m), inches=FALSE, add=T)
	# text(bSom$grid$pts[bins,] + rnorm(2*dim(x)[1],0,0.1), labels=rownames(data), cex=0.5, col=c("red","blue")[data$loc])
	points(bSom$grid$pts[bins,] + rnorm(2*dim(x)[1],0,0.1),...)
	text(rep(0, n), seq(1, n), labels = seq(1, n*m, by = n), cex = 1.5)
	text(rep(m+1, n), seq(1, n), labels = seq(n, n*m, by = n), cex = 1.5)
	list(Kohonen = bSom, classifier = cbind.data.frame(Ids = rownames(data), class = class))

}
