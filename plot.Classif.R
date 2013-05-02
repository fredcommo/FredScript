plot.Classif<-function(values, ord=TRUE, title=NULL, PCH = 19, CEX = 1.25, Ylab = "")

# Recherche des classes dans la série X et représente les mélanges de densité

# values 	: vecteur de valeurs
# infos 	: fichier info sur la variable à classifier
# ord 	: représentaion par valeurs croissantes

#	!!!! IMPORTANT !!!!
#     Nécessite library(mclust)

{
	x <- values
	n <- length(x)

	# Recherche d'un mélange gaussien
	BIC.x<-mclustBIC(x)
	x.model<-mclustModel(x,BIC.x)
	nclass<-length(x.model$parameters$mean)
	grps.x<-rep(1,length(x))

	# Classifieur et calcul de densité (cas où nclass>1)

	if (nclass>1)
	{
		x.class<-cbind.data.frame(x.model$z,class=0)

		for (i in 1:nrow(x.class))
			x.class$class[i]<-which.max(x.class[i,1:nclass])

		grps.x<-x.class$class
		table(grps.x)

		m<-x.model$parameters$mean
		s<-x.model$parameters$variance$sigmasq
		p<-x.model$parameters$pro
	
		X<-rnorm(n*p[1], mean=m[1], sd=sqrt(s[1]))

		for (i in 2:nclass)
			{
			if (length(s)==1)	X<-c(X, rnorm(n*p[i], mean=m[i], sd=sqrt(s)))
			else X<-c(X, rnorm(n*p[i], mean=m[i], sd=sqrt(s[i])))
			}
		dens<-density(X)
	}

	if (nclass==1) dens<-density(x)


	# Représentation du nuage de points et des densités
	n<-length(x)
	q05<-quantile(x, 0.05)	#*0.9
	q95<-quantile(x, 0.95)	#*1.1
	med <- quantile(x, 0.5)
	seq<-seq(0,1,by=1/(nclass+1))[2:(nclass+1)]
	seq1.x<-rnorm(length(which(grps.x==1)), mean=seq[1], sd=0.05)

	for (i in 2:nclass)
		seq1.x<-c(seq1.x,rnorm(length(which(grps.x==i)), mean=seq[i], sd=0.075))

	seq2.x<-seq(0,1,by=1/(n-1))

	grps.ord<-grps.x
	x.ord<-x

	if (nclass==1)
		{
			seq1.x<-seq2.x
			ord=FALSE
		}

	if (ord)
		{
			grps.ord<-grps.x[order(x)]
			x.ord<-sort(x)
		}

		# Représentation du nuage
		
		plot(x.ord~seq1.x, col=c("lightblue4", "orangered", "seagreen", "purple2", "burlywood2")[grps.ord],
			pch = PCH, cex = CEX, ylim = range(q05*0.9, q95*1.1), xlim=range(0,1), ylab= Ylab, xlab= "Fn(x)",
			main = paste("Values (per exp.)", title, sep=" : "))

		# Représentation des quantiles 0.05, 0.95 et médiane. Utilise les quantiles du dataset complet, donc à calculer avant !
		points(rep(q95, length(seq2.x))~seq2.x, cex=0.25, pch=25, col="grey")
		points(rep(q05, length(seq2.x))~seq2.x, pch=24, cex=0.25, col="grey")
		points(rep(med, length(seq2.x))~seq2.x, pch=22, cex=0.25)
		
		# Représentation de la courbe de densité
		lines(dens$x~I(dens$y/max(dens$y)), col="violetred4", lwd = 2)
		# lines(dens$x~dens$y)
		abline(v=0)
		
		# Legends
		legend("topleft", legend=paste("n =",n), cex=1, bty="n")
		legend(x=0.9, y=mean(q05), legend="q0.05", cex=0.75, text.col="black", bty="n")
		legend(x=0.9, y=mean(q95), legend="q0.95", cex=0.75, text.col="black", bty="n")
		legend(x=0.9, y=mean(med), legend="median", cex=0.75, text.col="black", bty="n")

	# Sortie groupes d'affectation
		
	grps.x
}
