permut3<-function(x,B)

# Version 3. Utilisation de sample(X)

# Permutations par groupe de 10 contre 10
{
log2<-x
B=B
tab.T3<-matrix(0,dim(log2)[1],B)
rownames(tab.T)<-rownames(log2)

fT3<-function(x)(t.test(x[1:10],x[11:20])$statistic)

	for (i in 1:B)
	{

		grp2<-log2[,sample(11 :100,10)]
		grp.test<-cbind(log2[,1:10],grp2)
		test<-sample(grp.test)

		tab.T3[,i]<-as.numeric(apply(test,1,fT3))
}
	T.bar3<-apply(tab.T3,1,mean)
	p.tab3<-dt(T.bar3,df=16)
	
	list(t.value=T.bar3,p.value=p.tab3)
}