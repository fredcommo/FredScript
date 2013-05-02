CEM.2<-function(x,k,nstart=10)

# x		: vecteur des valeurs
# k		: nombre de classes
# nstart	: nombre d'itérations

{
	  source("D:\\Stats\\Doc R\\Scripts R\\CEM.etape.E.2.R")

        ##################
        # initialisation #
        ##################

        #### Choix des centres initiaux

        x<-as.matrix(x);
        n<-dim(x)[1];
		
		# initialisation des centres       
		# MU<-matrix(0,k,1);
		# indexMU<-sample(1:n,k);
		# MU <- x[indexMU ,1 ];


	# Initialisation des centres .v2
	if (k==2) MU<-as.numeric(quantile(x,probs=c(0.2,0.8)))
	if (k==3) MU<-as.numeric(quantile(x,probs=c(0,0.5,1)))#c(min(x),median(x),max(x))
	if (k>3)
		{
			MU<-matrix(0,k,1);
			indexMU<-sample(1:n,k);
			MU <- x[indexMU ,1 ];
		}

        #initialisation des variances (tout le monde à 1)
        sigma<-rep(1,k);

        #initialisation des proportions
        prop<-rep(1/k,k);



        ##############
        # iterations #
        ##############

        for(i in 1:nstart)
        {
                	#### E

                	tik<-CEM.etape.E.2(x[,1],MU,sigma,prop)
			#tik<-tik+matrix(rnorm(n*k,sd=1e-5),n,k)

                	#### C
                	cluster<-apply(tik,1,which.max)       


                	#### M
                	groups<-split(x,cluster)

			MU<-lapply(groups,mean)
			MU<-unlist(MU) 
			sigma<-lapply(groups,var)
			sigma<-unlist(sigma)
                
        }

		prop<-table(cluster)/length(x) 
		prop<-as.numeric(prop)

        list(cluster=cluster,moyenne=MU,variance=sigma, proportion=prop)
}
