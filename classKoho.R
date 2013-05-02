classKoho<-function(x,y,n,m)
{
	lg<-length(x[,1]);
	pTEST<-floor(lg/6);
	pAPPR<-lg - pTEST;
	
	#Nombre de neurones
	nbn<-n*m;

	# Sélectionne des indiv. au hasard et extrait dans xmtest
	index.xm<-sample(1:lg,pAPPR)
	xmtest<-xm[index.xm,];
	
	# Extrait les classes norm/tum des patients sélectionnés 
	ytest<-as.matrix(y);
	ytest<- ytest[index.xm,];


	# Algorithme de Kohonen
	library(class);
	sg<-somgrid(n,m);
	radii<-c(seq(m,0,len=30),rep(0,10));
	batchSOM(xmtest,sg,radii)->alontest;

	# (knn1(alontest$codes,xmtest,1:nbn);)

	# Partitionne le groupe en fonction des affectations 
	bins<-as.numeric(knn1(alontest$codes,xmtest,1:nbn));

	# Construit le graphe
	plot(alontest$grid,type="n");
	symbols(alontest$grid$pts[,1],alontest$grid$pts[,2],circles=rep(0.4,nbn),inches=FALSE,add=TRUE);
	points(alontest$grid$pts[bins,]+rnorm(pAPPR*2,0,0.1),col=c("red","darkgreen")[y],pch=8,cex=2);

	### Construit la matrice des individus restants
	xmrest<-xm[-index.xm,];

	# Affecte les nouveaux patients et les place sur le graphe
	bins2<-as.numeric(knn1(alontest$codes,xmrest,1:nbn));
	points(alontest$grid$pts[bins2,]+rnorm(pTEST*2,0,0.1),col="blue",pch=8,cex=2);

	table(ytest,bins)->table_appr;
	vote_maj(table_appr);

}