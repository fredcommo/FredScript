CVknn<-function(x,kppv,V,...)	# funct(apprentissage,nb_ppv,nb_sousens.)
{
	library(class);
	rand<-sample(V,dim(x)[1],replace=T);
	knnres<-x[,1];
		# erreur<-rep(0,max(rand));
		# erreur_bar<-0;

	for(i in sort(unique(rand)))
		{
			test<-x[rand==i,];
			train<-x[rand!=i,];
			knnres[rand==i]<-knn(train[,-1],test[,-1],cl=train$letter,k=kppv);
			
		# Calcule l'erreur moyenne

		#	tab<-table(test$letter,knnres[rand==i]);
		#	nb<-sum(tab);			
		#	erreur[i]<-round(100*(sum(tab)-sum(diag(tab)))/sum(tab),2);
			
		}

		# erreur_bar<-1/max(rand)*sum(erreur);
		# cat("Erreur moyenne= ",erreur_bar,"%\n");
		# list(res=knnres,erreur_bar=erreur_bar);
	knnres;
}