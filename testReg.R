testReg<-function(X)
{

	#X<-rnorm(30,mean=100,sd=15);
	
	res<-matrix(0,2,1000);
	sig<-matrix(0,1,1000);

	for(i in 1:1000)
	{
		Y<- (1/6)*(X)+rnorm(30,mean=0,sd=2) -5;
		lm(Y~X)->g;
		g$coefficients->res[,i];
		(summary(g)$sigma)^2->sig[,i]
	}

	X1<-cbind(rep(1,30),X)
	vBh<-4*solve(t(X1)%*%X1) #variance attendue

	list(BetaMatrix=res,varBhat=vBh,varObs=cov(t(res)),sigma=sig);

}