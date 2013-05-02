na.remove<-function(x)
{
p<-dim(x)[2]
index<-rep(0,p)

	for (i in 1:p)
	{	tmp<-which(is.na(x[,i]))
		index<-c(index,tmp)
	}
	index<-sort(unique(index[index!=0]))
	x<-x[-index,]
	x
}
