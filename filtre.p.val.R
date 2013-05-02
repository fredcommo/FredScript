filtre.p.val<-function(x,seuil,freq)

#	renvoie les indices respectant les critères

#	x: data.frame des p.values
#	seuil: seuil de p.value à prendre en compte
#	freq: fréquence minimale respectant la p.value seuil.

{
	x<-as.data.frame(x)
	n<-dim(x)[1]
	p<-dim(x)[2]
	tmp<-x<=seuil
	freq.obs<-rep(0,n)
	
		for (i in 1:n)
			freq.obs[i]<-length(which(tmp[i,]))

	freq.obs<-freq.obs/p
	index<-which(freq.obs>=freq)
	index
}

B=100
n<-dim(x)[1]
p<-dim(x)[2]
res.t<-matrix(0,n,B)
#colnames(res)<-c("index","t.val","p.val")
#res[,1]<-index2
val.tmp<-rep(0,round(p*2/3))
#t.val<-rep(0,B)
for (i in 1:n)
{
	probe<-tmp2[i,]

B=1000
n<-dim(x)[1]
p<-dim(x)[2]
res.t<-matrix(0,n,B)
val.tmp<-rep(0,round(p*2/3))
	
		# effectue B répétitions du test
		for (j in 1:B)
		{ 	
			sample(1:p,round(p*2/3))->index.tmp
			val.tmp<-x[,index.tmp]
			means<-rowMeans(val.tmp)
			s<-sd(t(val.tmp))
			res.t[,j]<-means/(s*sqrt(1/length(index.tmp)))
		}
		p.val<-2*(1-pt(abs(res.t),length(index.tmp)-1))
		mean.t<-rowMeans(res.t)
		seuil.min<-min(mean.t)*0.25
		seuil.max<-max(mean.t)*0.25
		plot(mean.t,pch=8,cex=0.1)
		abline(h=c(seuil.min,seuil.max),lty=3,col="red")

	res[i,2]<-mean(t.val)
	res[i,3]<-t.test(t.val)$p.value
}

plot(density(t.val))



