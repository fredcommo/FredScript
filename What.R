What<-function(train,...)

{

	n<-dim(train)[1];				# nombre de lignes (individus)

	p<-dim(train[,-1])[2];				# nombre de v.a.

	table_class<-table(train$letter);		# nombre d'indiv. par classe letter
	
	k<-dim(table_class);				# nombre de classes

	prop<-table_class/n				# proportions (nk/n)


		# centrage des groupes: (xi-µk)

	splitdata<-split(train,train$letter);		# regroupe par classe

	split_cr<-splitdata

		for(i in 1:length(splitdata))
			split_cr[[i]]<-scale(splitdata[[i]][-1],scale=F);


		# Calcul des variances sigmak de chaque classe k: sigma_k=(1/nk)*t(xi-µk)(xi-µk)

	sigmak<-list(matrix(0,p,p));

		for(i in 1:k)
		sigmak[[i]]<-(1/table_class[i])*t(as.matrix(split_cr[[i]]))%*%(as.matrix(split_cr[[i]]));



		# Estimation de la variance intra-classe: Somme(k)[(nk/n)*sigma_i]

	W_hat<-matrix(0,p,p);
	
		for(i in 1:k)
		W_hat<-W_hat+prop[i]*sigmak[[i]];

	W_hat;


}