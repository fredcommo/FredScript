kinome.loess<-function(X, span=span, graph.title)
{
	x<-X$log1
	y<-X$log2
	#QC<-rep(NA,length(x))

	xbar<-mean(x,na.rm=TRUE)
	ybar<-mean(y,na.rm=TRUE)

	xi.squared<-(x-xbar)^2
	Sum.xi.squared<-sum(xi.squared,na.rm=TRUE)

	loess.res<-loess(I((x+y)/2)~x,span=span)
	true.mean<-(x+y)/2

	#new.x<-(x+y)/2
	#new.xbar<-mean(new.x)

	#new.xi.squared<-(new.x-new.xbar)^2
	xi.squared<-(x-xbar)^2

	#new.sum.squared<-sum(new.xi.squared)
	sum.squared<-sum(xi.squared)	

	#mean.fit<-predict(loess.res,newdata=new.x,se=TRUE)

	#na.index<-which(is.na(mean.fit$fit))
	#	if (length(na.index>0)) 
	#		{
	#		mean.fit$fit[na.index]<-loess.res$fitted[na.index]
	#		mean.fit$se[na.index]<-loess.res$residuals[na.index]
	#		}


	s2<-loess.res$s^2
	df<-length(x)-1
	t<-qt(0.025,df,lower.tail=F)

	#yi<-mean.fit$fit
	yi<-loess.res$fitted
	#yi.squared<-mean.fit$se
	yi.squared<-loess.res$residuals	

	CI<- rep(0,length(x))

	for (i in 1:length(x))
		#CI[i]<-t*s2*sqrt(1+1/length(!is.na(x))+(new.x[i]-new.xbar)^2/new.sum.squared)
		CI[i]<-t*s2*sqrt(1+1/length(!is.na(x))+(x[i]-xbar)^2/sum.squared)
		

	plot(y~x, pch=19, cex=0.5, col="blue",
		main=paste("Duplicates correlation", graph.title,sep="\n"), xlab="value1", ylab="value2")
	points(true.mean~x, pch=3, col="lightblue3", cex=0.4)
	lines(sort(yi)~sort(x),col="red",lty=1)
	lines(x=sort(x),y=sort(yi+CI), col="black")
	lines(x=sort(x),y=sort(yi-CI), col="black")

	

	#legend("topleft",legend=paste("Adj.R²",round(adj.r,4),sep=" = "),
	#					bty="n",cex=0.75,text.col="darkblue")
	
	#QC<-ifelse((fit.res$fit<val.min | y>val.max),"flag","OK")

	#QC[which(is.na(QC))]<-"suppr"
	#QC<-ifelse(x<=0 | y<=0,"suppr",QC)
	
	#list(fit=mean.fit$fit,resid=mean.fit$se)	# QC=QC
	list(fit=yi,resid=yi.squared)
}
