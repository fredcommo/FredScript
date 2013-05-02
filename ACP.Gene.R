
ACP.Gene <- function (eset, first=1, last=3, Pmax = 0.999){

# eset : data set with probes in rows, and exp in columns
# first, last : the first and last principal components to consider
# Pmax : the p-value to consider to define the Chisquared quantile.

	eset <- as.data.frame(eset)

# vérification
	prcomp(scale(eset))->acp1			# genes = obs ; exp = variables
	prcomp(scale(t(eset)))->acp2			# exp = obs ; genes = variables

	X <- as.data.frame(acp1$x[,first:last])
	X <-scale(X)
	X <- as.data.frame(X)

	Q <- qchisq(p = Pmax, df = ncol(X))
	D<-apply(X^2, 1, sum)
	pt.col<-ifelse(D>Q,"red","grey")

	inform<-which(pt.col=="red")

	sub1<-eset[inform,]
	acp3 <- prcomp(t(sub1))	# exp = obs ; genes = variables
	acp.mat <- acp3$x[,1:2]

	output <- list(Select = sub1, Inform = inform)

	par(mfrow=c(2,2))
	plot(acp1$x[,1:2], pch=19, cex=0.25, col="grey", asp=1, main="Probes")
	plot(acp2$x[,1:2], pch=19, cex=1.2, main="Initial PCA")
	plot(acp1$x[,1:2], cex=0.25, col="blue", pch=19, asp=1, main=paste("Informative probes;", Pmax))
	points(acp1$x[,2]~acp1$x[,1], cex=0.2, pch=19, col=pt.col)

	plot(acp.mat, pch=19, cex=1.2, main="PCA with informative probes")

	output
}