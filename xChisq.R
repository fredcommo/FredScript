xChisq<-function(X)
{
	p<-dim(X)[2]
	p.val<-matrix(0,p,p)
	for (i in 1:p)
	{ 
		for (j in 1:p)
		{
			p.val[i,j]<-chisq.test(X[,i],X[,j],
				rescale.p=T,simulate.p.value=T)$p.value
		}
	}
	p.val<-signif(p.val,4)
	p.val<-round(p.val,4)
	rownames(p.val)<-colnames(X)
	colnames(p.val)<-colnames(X)
	p.val<-as.data.frame(p.val)
	p.val
}