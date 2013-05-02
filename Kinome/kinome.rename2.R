kinome.rename2<-function(X,File.infos)
{
	
	# create a data frame to identify wells
	# data.frame info sequences : fichier de sortie

	result<-cbind.data.frame(File.infos[,2:5],value1=0,value2=0)
			
	
	n<-dim(result)[1]

	for (i in 1:n)
	{
		well<-result$well[i]
		tmp<-X[which(X$well==well),]
		result$value1[i]<-as.numeric(tmp$value[1])
		result$value2[i]<-as.numeric(tmp$value[2])
	}
		
	result
}