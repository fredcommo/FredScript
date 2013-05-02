kinome.means<-function(X)
{
	p<-dim(X)[2]
	lev<-levels(X$well)
	n.lev<-length(lev)-1
	
# create a data frame
	
	resp<-X[1,]				# data.frame info séquences

	i=1
	
	for (j in 1:16)
		for (k in 1:24)
		{
			resp[i,1]<-paste(LETTERS[j],k,sep="")
			i=i+1
		}	


	for (i in 1:n.lev)
	{
		well<-resp$well[i]
		resp$value[i]<-mean(X[which(X$well==well),]$value)
		resp[i,-c(1:2)]<-X[which(X$well==well),3:p][1,]
	}

	resp
}