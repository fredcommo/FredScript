# bimodality index : utiliser les kmeans (c=2) ou mclust

# in case of failure, try(kmeans()) return an error message
# this may be due to absence of distinct points : all NA,...

BIndex <- function(x){

	x <- as.numeric(x)
	test <- try(kmeans(x, centers = 2), silent = T)

	if(class(test)!="try-error"){
		km1 <- kmeans(x, centers = 2)
		clust <- km1$cluster
		mu <- km1$centers
		ni <- km1$size
		n <- sum(ni)

		S1 <- var(x[clust==1], na.rm = T)
		S2 <- var(x[clust==2], na.rm = T)

		s <- ((ni[1]-1)*S1 + (ni[2]-1)*S2)/(n-2)
		p <- ni[1]/n
		delta <- abs(mu[1] - mu[2])/sqrt(s)
		BI <- ((p*(1-p))^(1/2))*delta		# suggested cutoff : BI>=1.1 
		}
	else BI <- NA
	return(BI)	
	}
