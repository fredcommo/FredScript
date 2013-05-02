kinome.fitlm<-function(X)

# X	: $Values from kinome.QC or normalization

{
	result<-cbind.data.frame(X$Values,fit.val=0,fit.sd=0)	

	x<-X$Values$value1
	y<-X$Values$value2

	lm.val<-lm(y~x)
	result$fit.val<-predict(lm.val)

	xbar<-mean(x)
	ybar<-mean(y)

	adj.r<-summary(lm.val)$adj.r.squared
	residuals<-lm.val$residuals
	S.res<- summary(lm.val)$sigma
	df<-lm.val$df.residual
	
	result$fit.sd<-S.res*sqrt(1+1/length(x)+(residuals-xbar)^2/sum((residuals-xbar)^2))

	result

}