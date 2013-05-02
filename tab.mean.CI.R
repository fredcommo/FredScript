tab.mean.CI<-function(X,crit)
{
	crit<-as.data.frame(crit)
	ligne=0
	
	# compter les critères dans crit: =c
	c<-dim(crit)[2]
	
	for (i in 1:c){
	ligne=ligne+length(levels(crit[,i]))
	}
	
	# compter les variables dans X
	p<-dim(X)[2]
	col<-p+2

	res<-matrix(0,ligne,col)
	colnames(res)[3:(col)]<-colnames(X)
	colnames(res)[1:2]<-c("Critère","niveau")
	
	# sélectionner le critère
	index=1
	for (j in 1:c)
	{
		# compter les niveaux du critère: =l
		l<-length(levels(crit[,j]))
	
		# pour chaque niveau:
		for (k in 1:l)
		{
			tmp<-X[crit[,j]==levels(crit[,j])[k],]
			mu<- signif(mean(tmp,na.rm=TRUE),4)
			sd<-sd(tmp,na.rm=TRUE)
			n<-colSums(!is.na(tmp))
			t<-qt(0.975,n)
			ci<- signif(t*sd*sqrt(1/n),3)
			res[index,3:(p+2)]=paste(mu,ci,sep=" +/-")
			res[index,1]<-colnames(crit)[j]
			res[index,2]<-levels(crit[,j])[k]
			index=index+1
		}
	}
	res<-as.data.frame(res)
	res
}
