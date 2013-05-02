permut.moy<-function(grpA,grpB,B)
{
n<-dim(grpA)[2]
m<-dim(grpB)[2]

index<-matrix(0,B,dim(grpA)[1])
stat<-matrix(0,B,dim(grpA)[1])
colnames(stat)<-rownames(grpA)
p.val<-rep(0,dim(grpA)[1])
p.adj<-rep(0,dim(grpA)[1])

class<-c(rep(0,n),rep(1,n))

	for (i in 1:B)
	{
		grp2<-grpB[,sample(1:m,n)]
		grp.test<-cbind(grpA,grp2)
		res<-mt.maxT(grp.test,class,B=1000)
		res<-res[order(res$index),]
		stat[i,]<-res$teststat
	}
Stat<-apply(stat,2,mean)

p.val<-dt(Stat,df=16)
p.adj<-p.adjust(p.val,method="fdr")

cbind(stat=Stat,rawp=p.val,adjp=p.adj)->tableau
as.data.frame(tableau)->tableau
tableau
}



