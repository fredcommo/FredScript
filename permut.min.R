# Version 7. Utilisation de mt.maxT en boucle avec recherche de min(adjp)

permut.min<-function(grpA,grpB,B)
{
n<-dim(grpA)[2]
m<-dim(grpB)[2]

rawp<-matrix(0,B,dim(grpA)[1])
colnames(rawp)<-rownames(grpA)
p.adj<-matrix(0,B,dim(grpA)[1])
stat<-matrix(0,B,dim(grpA)[1])



class<-c(rep(0,n),rep(1,n))

	for (i in 1:B)
	{
		grp2<-grpB[,sample(1:m,n)]
		grp.test<-cbind(grpA,grp2)
		res<-mt.maxT(grp.test,class,B=1000)
		res<-res[order(res$index),]
		rawp[i,]<-res$rawp
		p.adj[i,]<-res$adjp
	}
	Rawp<-apply(rawp,2,min)
	Adj<-apply(p.adj,2,min)

	cbind(rawp=Rawp,adjp=Adj)->tableau
	as.data.frame(tableau)->tableau
	tableau
}
