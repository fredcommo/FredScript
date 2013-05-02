Bhat<-function(train,...)

{

	n<-dim(train)[1];				# nombre de lignes (individus)

	p<-dim(train[,-1])[2];				# nombre de paramètres

	table_class<-table(train$letter);		# nombre d'indiv. par classe letter
	
	k<-dim(table_class);				# nombre de classes

	prop<-table_class/n;				# proportions (nk/N)


		# Calcule le vecteur moyennes 

	mean(train[,-1])->MU


		# Calcul du vecteur moyennes pour chaque classe k

	mu_class<-matrix(0,k,p);
	mu_class<-aggregate(train[,-1],list(letter=train$letter),mean);

		
		# Calcul de la variance inter_classe : Somme(k)[(nk/n)*t(µk-µ)(µk-µ)]

	B_hatk<-list(matrix(0,p,p));

		for(i in 1:p)
		B_hatk[[i]]<-prop[i]*t(as.matrix(mu_class[i,-1]-MU))%*%as.matrix(mu_class[i,-1]-MU);

	B_hat<-matrix(0,p,p);
		
		for(i in 1:p)
		B_hat<-B_hat+B_hatk[[i]];

	B_hat;

}

		