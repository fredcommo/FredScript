CEM.etape.E.2<-function(x,MU,sigma,prop)
{
        n<-length(x)
        k<-length(MU)
        Xi<-matrix(x,n,k)
        sigmak<-t(matrix(sigma,k,n))
        pk<-t(matrix(prop,k,n))
        MUk<-t(matrix(MU,k,n))

        logtik1<- log(pk) + log(1/sqrt(2*pi*sigmak^2)) - (1/(2*sigmak^2))*(Xi-MUk)^2
	  logtik2<-log(apply(exp(logtik1),1,sum))

        logtik<- logtik1 - logtik2
        
        return(exp(logtik))
}
