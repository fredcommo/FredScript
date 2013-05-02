kinome.sam<-function(X,Y)
{
	p<-dim(X)[2]
	lev<-levels(X$well)
	n.lev<-length(lev)-1
	

	result<-cbind.data.frame(well=X[1,1],t=0,p.value=0,adj.p=0,X[1,3:p])
	#info<-X[1,3:p]

	s0.X<-sd(X[which(X$well=="Blank"),]$value)
	s0.Y<-sd(Y[which(Y$well=="Blank"),]$value)

	i=1
	
	for (j in 1:16)
		for (k in 1:24)
		{
			result[i,1]<-paste(LETTERS[j],k,sep="")
			i=i+1
		}	


	for (i in 1:n.lev)
	{
		well<-result$well[i]

		mu1<-mean(X[which(X$well==well),]$value)
		sd1<-sd(X[which(X$well==well),]$value)

		mu2<-mean(Y[which(Y$well==well),]$value)
		sd2<-sd(Y[which(Y$well==well),]$value)
		
	# A vérifier
		t<-(mu1-mu2)/(sqrt(1/2*(sd1^2+s0.X^2+sd2^2+s0.Y^2)))
		p.value<-1-pt(abs(t),df=n.lev-2)

		result$t[i]<-signif(t,4)
		result$p.value[i]<-signif(p.value,3)
		result[i,5:(p+2)]<-X[which(X$well==well),3:p][1,]
	}
	
	adj.p<-p.adjust(result$p.value,method="BH")
	result$adj.p<-signif(adj.p,3)
	
	result
}
