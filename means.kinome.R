kinome.means<-function(X)
{

	lev<-levels(X[,1])
	n.lev<-length(lev)-1
	
# create a data frame

	resp<-matrix(0,n.lev,2)
	resp<-as.data.frame(resp)
		i=1
		for (j in 1:16)
			for (k in 1:24){
				resp[i,1]<-paste(LETTERS[j],k,sep="")
				i=i+1
			}	
	Blank<-mean(X[which(X[,1]=="Blank"),2])

	for (i in 1:n.lev)
	{
		well<-resp[i,1]
		resp[i,2]<-mean(X[which(X[,1]==well),2])-Blank
	}
	colnames(resp)<-c("well","value")
	resp
}