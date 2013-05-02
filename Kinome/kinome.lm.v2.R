kinome.lm.v2<-function(X)
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
		t<-qt(0.025,df,lower.tail=F)

	yi<-as.numeric(predict(lm.res, na.action = na.exclude))
	yi.squared<-(yi-ybar)^2

	CI<- rep(0,length(x))

	for (i in 1:length(x))
		CI[i]<-t*S.residuals*sqrt(1+1/length(!is.na(x))+(x[i]-xbar)^2/Sum.xi.squared)


	val.max<-yi+CI
	val.min<-yi-CI
	QC<-ifelse((y<val.min | y>val.max),"flag","OK")

	plot(y~x,pch=20,cex=0.75,col="blue",
		main="Duplicates correlation",xlab="value1",ylab="value2")
	abline(lm.res,col="red",lty=2)
	lines(x=x,y=c(yi+CI), col="grey")
	lines((yi-CI)~x, col="grey")
	legend("topleft",legend=paste("Adj.R²",round(adj.r,4),sep=" = "),
						bty="n",cex=0.75,text.col="darkblue")
	QC[which(is.na(QC))]<-"suppr"
	QC<-ifelse(x<=0 | y<=0,"suppr",QC)
	return (QC)

}
