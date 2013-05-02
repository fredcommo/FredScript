kinome.BlankCor<-function(X)
{
	result<-X
	n<-dim(result)[1]	
	Blank<-mean(X[which(X$kinase=="neg_control"),]$value)
		
	for (i in 1:n)
		result$value[i]<-result$value[i]-Blank
	
	result
}

