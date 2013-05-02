classKoho<-function(x,y,n,m)
{
	lg<-length(x[,1]);
	pTEST<-floor(lg/6);	#Ensemble de test = 1/6
	pAPPR<-lg - pTEST;

	
	#Nombre de neurones entrÈs par l'utilisateur
	nbn<-n*m;

	# Sélectionne des indiv. au hasard et extrait dans x_appr
	index.x<-sample(1:lg,pAPPR)
	x_appr<-x[index.x,];

	
	
	# Extrait les classes normal/tumeur des patients sélectionnés 
	y_appr<-as.matrix(y);
	y_appr<- y_appr[index.x,];


	# Algorithme de Kohonen sur les patients pour l'apprentissage
	library(class);
	sg<-somgrid(n,m);
	radii<-c(seq(m,0,len=30),rep(0,10));
	batchSOM(x_appr,sg,radii)->alontest;

	# Partitionne le groupe en fonction des affectations 
	bins<-as.numeric(knn1(alontest$codes,x_appr,1:nbn));

	# Construit le graphe
	plot(alontest$grid,type="n",main="Répartition des individus sur la carte");
	symbols(alontest$grid$pts[,1],alontest$grid$pts[,2],circles=rep(0.5,nbn),inches=FALSE,add=TRUE);
	points(alontest$grid$pts[bins,]+rnorm(pAPPR*2,0,0.1),col=c("red","darkgreen")[y],pch=8,cex=2);

	### Construit la matrice des individus restants
	xrest<-x[-index.x,];
	yrest<-as.matrix(y[-index.x]);

	# Affecte les nouveaux patients et les place sur le graphe
	bins2<-as.numeric(knn1(alontest$codes,xrest,1:nbn));
	points(alontest$grid$pts[bins2,]+rnorm(pTEST*2,0,0.1),col="blue",pch=8,cex=2);

	# Résultats
		# Table des neurones
		table(y_appr,bins)->table_appr;

		# Table de classement des patients testÈs
		table(rownames(xrest),bins2)->table_rest;

		# Table des vraies classes pour les patients testÈs
		cbind(rownames(xrest),yrest)->vraiclasse;
	
						#vote_maj(table_appr)->classes;

	# Affectation d'une classe à chaque neurone par vote majoritaire

		#res<-matrix("Indet",length(table_appr[1,]),1);	###  Indet par dÈfaut
		res<-matrix("Indet",nbn,1);
		#rownames(res)<-rep("",length(table_appr[1,]));
		rownames(res)<-rep("",nbn);

			for(i in 1:length(table_appr[1,]))
			{
					# Renomme les lignes par le N° du neurone
					rownames(res)[as.numeric(colnames(table_appr))[i]]<-colnames(table_appr)[i]
			
					if(table_appr[1,i]<table_appr[2,i])
						res[as.numeric(colnames(table_appr))[i],]<-rownames(table_appr)[2];

					if(table_appr[1,i]>table_appr[2,i])
						res[as.numeric(colnames(table_appr))[i],]<-rownames(table_appr)[1];
			}
	
	


	# Construit la matrice des résultats Classement vs rÈel.

		classif<-matrix("",length(xrest[,1]),3);
		cbind(rownames(xrest),bins2)->classif;	

		# Remplace le numéro de neurone par sa classe
			classif[,2]<-res[bins2,];
			cbind(classif,vraiclasse[,2])->classif;
			colnames(classif)<-c("Patient","classé","réel");


	# Calcule le % d'erreur de classement
		prop<-0;
		propindet<-0;
		for(i in 1:length(xrest[,1]))
			{
			if(classif[i,2]!=classif[i,3] && classif[i,2]!="Indet")
				prop<-prop+1;
			if(classif[i,2]=="Indet")
				propindet<-propindet+1;
			}
			prop<-(prop/length(xrest[,1]))*100;
			propindet<-(propindet/length(xrest[,1]))*100;

	# Sorties
	list(modele=alontest, Affectation=bins, table_appr=table_appr, Classe.neurones=res, classement=table_rest,Statut.reel=vraiclasse,classification=classif, proportion.erreur=prop, proportion.indet=propindet);
}