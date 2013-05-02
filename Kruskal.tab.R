Kruskal.tab<-function(class,val)

# Renvoie un tableau de p-values d'un test de rangs (Kruskal-wallis).

# 	class: vecteur de type "facteur" contenant les modalités à discriminer
#	val: data.frame contenant les variables.
{
	Res<-data.frame(name=colnames(val),df=0,Chisq=0,p.val=0)

	for (i in 1:dim(val)[2])
	{
		k<-kruskal.test(val[,i]~class)
		Res[i,]<-data.frame(name=colnames(val)[i],df=k$parameter,
		Chisq=signif(k$statistic,3),p.val=signif(k$p.value,3))
	}
	Res
}
