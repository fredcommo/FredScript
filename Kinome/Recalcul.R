Recalcul<-function(X)

# mu.Log ans var.log recalculation after normalization.
# The function also redefine the QC status

{
	new<-cbind.data.frame(X[,-9],mu.log=0,var.log=0,QC=X[,9])

	new$mean<-apply(new[,5:6],1,mean)
	new$sdev<- apply(new[,5:6],1,sd)

	n<-dim(X)[1]

	mu<- new$mean
	s2<- new$sdev^2

	new$var.log<-log10(s2/mu^2+1)

	for (i in 1:n)
			{
			if (new$value1[i]<0 | new$value2[i]<0) 
				{
					new$mean[i]<-mu[i]<-1
					new$QC[i]<-"flag"
				}
			}



	for (i in 1:n)
	{
		if (new$mean[i]>1)
		{
		new$mu.log[i]<-log10(mu[i])-1/2*log10(s2[i]/mu[i]^2+1)
			}			
		else
		{
			new$mu.log[i]<-0
			}
	}
	new
}

