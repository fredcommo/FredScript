Recalcul.2<-function(X)

# mu.Log ans var.log calculation after normalization.
# The function also redefine the QC status

{
	new<-cbind.data.frame(X[,-9],mu.log=0,var.log=0,QC=X$QC)

	n<-dim(X)[1]
	mu<- new$mean
	s2<- new$sdev^2
	new$var.log<-log10(s2/mu^2+1)

	for (i in 1:n)
		if (new$QC[i]=="suppr" | new$QC[i]=="flag") new$mean[i]<-mu[i]<-1

	new$mu.log<-ifelse(new$mean>1,log10(mu)-1/2*log10(s2/mu^2+1),0)

	new
}

