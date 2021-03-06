tab.mean.sd<-function(X,crit)
{
	crit<-as.data.frame(crit)
	max=0
	
	# compter les crit�res dans crit: =c
	c<-dim(crit)[2]
	
	# compter les variables dans X
	p<-dim(X)[2]

	for (i in 1:c){
	max=max+length(levels(crit[,i]))
	}
	res<-matrix(0,max,p+2)
	colnames(res)[3:(p+2)]<-colnames(X)
	colnames(res)[1:2]<-c("Crit�re","niveau")
	
	# s�lectionner le crit�re
	index=1
	for (j in 1:c)
	{
		# compter les niveaux du crit�re: =l
		l<-length(levels(crit[,j]))
	
		# pour chaque niveau:
		for (k in 1:l)
		{
			tmp<-X[crit[,j]==levels(crit[,j])[k],]
			mu<- signif(mean(tmp,na.rm=TRUE),4)		
			sd<- signif(sd(tmp,na.rm=TRUE),3)

			res[index,3:(p+2)]=paste(mu,sd,sep=" sd=")
			res[index,1]<-colnames(crit)[j]
			res[index,2]<-levels(crit[,j])[k]
			index=index+1
		}
	}
	res<-as.data.frame(res)
	res
}
