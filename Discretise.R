discretise<-function(x,nclass)

# discrétise l'ensemble des valeurs en n classes d'effectifs égaux

{
	seq<-function(dimen,nclass)
	{	
	res<-rep(0,(nclass-1))
	for (i in 1:(nclass-1))
	res[i]<-round(length(sort(tmp))/nclass)*i
	res
	}
	
	p<-dim(x)[2]
	
	for (i in 1:p)
		{
		tmp<-x[,p]
		dimen<-length(sort(tmp))
		sequence<-seq(dimen,nclass)
		c=1
		class<-rep(0,(nclass-1))
		for (j in sequence)
			{
			class[c]<-sort(tmp)[j]
			c=c+1
			}
		
		tmp[tmp<class[1]]<-paste("lev",1,sep="")
		for (k in 1:length(class))
		{
		tmp[tmp<class[k]]<-paste("lev",k,sep="")
		}
		x[,p]<-tmp
		}
		x
}