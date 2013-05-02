kinome.lm.v3<-function(X)
{
	x<-X$value1
	y<-X$value2
	QC<-rep(NA,length(x))

	xbar<-mean(x,na.rm=TRUE)
	ybar<-mean(y,na.rm=TRUE)

	xi.squared<-(x-xbar)^2
	Sum.xi.squared<-sum(xi.squared,na.rm=TRUE)

	lm.res<-lm(y~x,na.action = na.exclude)
		adj.r<-summary(lm.res)$adj.r.squared
		S.residuals<- summary(lm.res)$sigma
		df<-lm.res$df.residual
		# t<-qt(0.025,df,lower.tail=F)
	
	ratio<-x-y		# x,y are in log2

	mu<-mean(ratio)
	s<-sd(ratio)
	n<-length(ratio)

	QC<-rep(NA,n)

	# outliers threshold are defined as  mu +/- 2*sigma
	#val.max<-mu+2*s
	#val.min<-mu-2*s		

	# outliers threshold are defined as:h=1.5(q3-q1), low=q1-h; high=q3+h
	q1<-quantile(ratio,probs=0.25)
	q3<-quantile(ratio,probs=0.75)	
	h<-1.5*(q3-q1)
	val.max<-q3+h
	val.min<-q1-h

	QC<-ifelse((ratio<val.min | ratio>val.max),"flag","OK")
	QC<-ifelse(x<=0 | y<=0,"suppr",QC)
	
	plot(y~x,pch=20,cex=0.75,col="blue",
		main="Duplicates correlation",xlab="value1",ylab="value2")
	abline(lm.res,col="red",lty=2)
	legend("topleft",legend=paste("Adj.R²",round(adj.r,4),sep=" = "),
						bty="n",cex=0.75,text.col="darkblue")
		
	return (QC)

}
