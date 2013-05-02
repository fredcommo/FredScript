kinome.normalize<-function(X)

# X : concaténation des duplicates de chaque fichier

{
	library(affy)

	n<-dim(X)[2]/2
	p<-dim(X)[1]

	all1<-X[,1]
	all2<-X[,2]
	for (i in 2:n)
		{
		all1<-cbind(all1,X[,2*i-1])
		all2<-cbind(all2,X[,2*i])
		}
	all<-rbind(all1,all2)

	norm<-normalize.quantiles(all)

	for (i in 1:n)
		{
		X[,2*i-1]<-norm[(1:p),i]
		X[,2*i]<-norm[(p+1):(2*p),i]
		}
	X
}