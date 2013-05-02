# Function start

KM.Dyn <- function (X, Niter = 10, Plot = TRUE) {

# X : A matrix or data.frame, of dim = n,2, to clusterize
# Niter : number of iterations.
# Plot : should the graph be plotted at each iteration?

	N <- nrow(X)
	Kmax <- floor(sqrt(N))

	# plot(X, main="Analysis in progress...")

	TIC <- matrix(0, Niter, Kmax)
	colnames(TIC) <- paste("ClustBy", seq(1, Kmax), sep="")

	for (iter in 1:Niter){
		clust <- rep(1, N)
		size <- N

		for(k in 2:Kmax){
			KM <- kmeans(X, k, nstart=100)
			new.clust <- KM$cluster

			if (Plot){
				plot(X, col=c(1:k)[KM$cluster], pch=c(1:k)[KM$cluster], cex=0.8, main = "In progress...")
				title(sub=paste("Iteration : ", iter, ",","k =",k))
				}

			prop.tab <- table(new.clust, clust)/colSums(table(new.clust, clust)); prop.tab		

			# D = -2 si perdu, 2 sinon
				D <- rep(1, k-1)
				for(i in 1:k){
					# D[i] <- length(unique(clust[which(new.clust==i)]))
					if(length(which(prop.tab[i,]!=0))>1) D[which(prop.tab[i,]!=0)] <- (-1)
				}
				# D

			# pjk = prop des élements hérités par cluster le jème cluster issu du cluster i parent.

			M <- rep(0, k-1)
			for (i in 1:length(M)){
				p <- prop.tab[which(prop.tab[,i]!=0),i]
				M[i] <- sum(-p*log(p))
				}
				# M
	
			g = D*M ; g

			Nk <- size
	
			NIFTI = sum(Nk/N*g); NIFTI

			TIC[iter, k] = TIC[iter, k-1] + NIFTI; TIC

			clust <- new.clust
			size <- KM$size

		}
	}

	M.TIC <- mean(as.data.frame(TIC)); M.TIC

	best.k <- which.max(M.TIC); cat("Optimal clustering =", best.k, "groups","\n")

	best.KM <- kmeans(X, best.k, nstart=100)
	best.clust <- best.KM$cluster

	plot(X, col=c(1:best.k)[best.clust], pch=c(1:best.k)[best.clust], cex=0.8)
	title (main = paste("Optimal clustering =", best.k, "groups"))

	list(M.TIC = M.TIC, Best.k = best.k, Cluster = best.clust)
}

# Function end

# Values:
# M.TUC : mean of TICs for k=1 to Kmax
# Best.k : the best number of clusters
# Cluster : the cluster Id for each observation.
