
boxplot.colors<-function(X, g=1, a=0.5) {

# X : values in any format.
# g, a : color (rainbow) parameters in [0,1]

	if(!is.data.frame(X)) X <- as.data.frame(X)

	X.med <- apply(X, 2, median, na.rm = TRUE)
	col.palette = rainbow(length(X.med), start=0.6, end=0.16, gamma=g, alpha=a)

	box.col= col.palette[as.numeric(rank(X.med))]
	boxplot(X, col = box.col, border = box.col, names = NULL,
			outcol = box.col, outpch = 19, outcex = 0.25)
}

