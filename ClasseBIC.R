ClasseBIC<-function(x)
{
        n<-dim(x)[1];
        p<-dim(x)[2];
        BIC1<-rep(0,n-1);
        BIC2<-rep(0,n-1);

	  source("D:\\SY09\\TP3\\kmeans_nstart.R");
        
	
        for(k in 2:(n-1))	# k =	nb de classes possibles

        {
                res<-kmeans_nstart(x,k,nstart=10);
                trS<-(1/n)*sum(res$withinss);
                lc<-(-1/2)*(n*p+n*p*(log(trS/p))) - n*log(k) + ((n*p)/2) * log(2*pi);
                BIC1[k]<- -2*lc 
                BIC2[k]<- (p*k+1)*(log(n));
        }
        list(BIC1,BIC2);
}