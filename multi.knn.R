multi.knn<-function(x,k,class=1,replic)

# Extrait un ens. appr et test n.fois (rep) et effectue knn avec k=k
# Le nombre d'erreur est calculé à chaque rep et ajouté dans un vecteur(err)


{	library(class)
	
	err<-rep(0,replic)

	for (i in 1:replic)
	{
		n<-dim(x)[1]
		index<-sample(1:n,round(n/3))

		train<-x[index,]
		test<-x[-index,]

		knn.res<- knn(train[,-class],test[,-class],cl=train[,class],k=k)
		tab<-table(test$Loc,knn.res)
	
		err[i]<-round(100*(sum(tab)-sum(diag(tab)))/sum(tab),2)
	}
	err

}