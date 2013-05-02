Effet<-function(data,perm)
{
	n<-dim(data)[1]
	names<-rownames(data)
	Effect<-rep(0,n)
	
	for (i in 1:n)
		{
		Eff<-rep(0,perm)
	
			for (j in 1:perm)
			{
				index1<-sample(1:10,5)
				index2<-sample(11:16,5)
				MU<-(1/10)*sum(data[i,c(index1,index2)])
				mu1<-(1/5)*sum(data[i,index1])
				Eff[j]<-mu1-MU
			}
		Effect[i]<-mean(Eff)
		}
	Effect<-data.frame(Effect)
	rownames(Effect)<-rownames(data)
	Effect
}