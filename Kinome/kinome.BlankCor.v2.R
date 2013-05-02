kinome.BlankCor.v2<-function(X,neg.ctrl="neg._control")
{

	result<-X
	n<-dim(result)[1]
	cor.values<-X[which(X$kinase==neg.ctrl),]$value

	m<-mean(cor.values)
	l<-length(cor.values)
	s<-sd(cor.values)
	t<-qt(0.975,l-1)
	IC<-t*s*sqrt(1/l)
	
	Blank<-m+IC
		
	for (i in 1:n)
		result$value[i]<-result$value[i]-Blank
	
	result$well<-as.factor(result$well)
	result$value[which(result$value<=0)]<-NA
	
	result
}

