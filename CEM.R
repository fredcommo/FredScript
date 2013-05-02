CEM<-function(x,k,nstart)
{
	#initialisation
	#choix des centres initiaux
	n<-dim(x)[1];
	p<-dim(x)[2];
	MU<-matrix(0,k,p);
	indexMU<-sample(1:n,k);
	MU <- x[indexMU , ];
}