CEM.3<-function(x,k,nstart=10)

# x		: dataset (ACP: PC1 to PC3)
# k		: nombre de classes
# nstart	: nombre d'itérations

# dimension = 3, remplacer par d dimensions

{
	x<-as.matrix(x)
	n<-dim(x)[1]
		
	# initialisation des centres       
	MU<-matrix(0,k,3)
	indexMU<-sample(1:n,k)
	MU <- x[indexMU , 1:3]
	

	#initialisation des variances (tout le monde à 1)
	sigma<-matrix(1,k,3)

	#initialisation des proportions
	prop<-rep(1/k,k)

for(i in 1:nstart)
        {
                	#### E step

                	{
        		n<-nrow(x)
        		k<-nrow(MU)
        		Xi<-matrix(x,n,ncol(x)*k)
        		
			sigma<-as.vector(sigma)
			sigmak<-matrix(sigma,n,ncol(x)*k, byrow=T)

			pk<-matrix(prop,n,ncol(x)*k)

			MU<-as.vector(t(MU))
 			MUk<-matrix(MU,n,ncol(x)*k, byrow=T)
				
        		logtik1<- log(pk) + log(1/sqrt(2*pi*sigmak^2)) - (1/(2*sigmak^2))*(Xi-MUk)^2
			logtik2<-log(apply(exp(logtik1),1,sum))

        		logtik<- logtik1 - logtik2
        
        		tik<-exp(logtik)
			}


                	#### C step
                	cluster<-apply(tik,1,which.max)       


                	#### M step
                	groups<-split(as.data.frame(x),cluster)

			MU<-lapply(groups,mean)
			MU<-t(matrix(unlist(MU),k,3)) 
			sigma<-lapply(groups,var)
			sigma<-t(matrix(unlist(sigma),k,3))
                

			prop<-table(cluster)/nrow(x) 
			prop<-as.numeric(prop)
			k<-length(prop)
		}
	graph3D.3(acp4, class1=as.factor(cluster), size=2)
	
	list(cluster=cluster,moyenne=MU,variance=sigma, proportion=prop)
}
