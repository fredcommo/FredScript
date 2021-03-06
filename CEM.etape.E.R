CEM.etape.E.2<-function(x,MU,sigma,prop)
{
        n<-length(x);
        k<-length(MU);
        Xi<-matrix(x,n,k);
        sigmak<-t(matrix(sigma,k,n));
        pk<-t(matrix(prop,k,n));
        MUk<-t(matrix(MU,k,n))

        logtik<- log(pk) + log(1/sqrt(2*pi*sigmak^2)) - (1/(2*sigmak^2))*(Xi-MUk)^2;

        logtik<- logtik - log(apply(exp(logtik),1,sum));
        
        return(exp(logtik));
}
