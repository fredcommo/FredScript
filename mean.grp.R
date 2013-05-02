mean.grp<-function(x,crit1,crit2,crit3,cut,...)

# calcule les moyennes par critère

{
l1<-length(levels(crit1))
l2<-length(levels(crit2))
l3<-length(levels(crit3))
p<-dim(x)[2]-3

Res<-matrix(0,l1*l2*l3,p)
n=1

	for (i in 1:l1)
		for (j in 1:l2)
			for (k in 1:l3)
			{ 	tmp<-x[x[,1]==levels(crit1)[i] 
					& x[,2]==levels(crit2)[j] 
					& x[,3]==levels(crit3)[k],-c(1:3)]
				if (nrow(tmp)>=cut) Res[n,]<-round(mean(tmp,trim=0.1,na.rm=TRUE),2)
				n=n+1
			}
	col1<-rep(levels(crit1),each=l2*l3)
	col2<-rep(levels(crit2),each=l3)
	col3<-rep(levels(crit3),l1*l2)
	Res<-cbind(col1,col2,col3,as.data.frame(Res))
	colnames(Res)<-colnames(x)
	index<-which(rowSums(Res[,-c(1:3)])==0)
	Res<-Res[-index,]
	heatmap(as.matrix(Res[,-c(1:3)]),col=cm.colors(100,1),Rowv=NA,
		scale="column",labRow=paste(Res[,1],Res[,2],Res[,3],sep="-"))
	Res
}
#Rowv=NULL,Colv=NULL

