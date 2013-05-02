formate<-function(X, span=0.75, graph.title=NULL)
{
	if (!is.factor(X$Flags))
	{
		X$Flags<-factor(X$Flags)
		levels(X$Flags)<-c("Flag","OK")
	}
	slide.lev<-levels(X$Slide)

	output<-cbind.data.frame(X[1,1:7], Intensity1=0, Intensity2=0, Mean=0, SDev=0,
					log1=0, log2=0, log.mean=0, mean.log=0, fit=0, resid=0)
	
	well.Id<-paste(LETTERS[1],seq(1,24),sep="")
	for (i in 2:16)
		well.Id<-c(well.Id,paste(LETTERS[i],seq(1,24),sep=""))
	
	n<-length(well.Id)
	
	k=1

	for (i in 1:length(slide.lev))
	{
		tmp<-X[which(X$Slide==slide.lev[i]),]

		for (j in 1:n)
		{
			tmp2<-tmp[which(tmp$Spot.Coord==well.Id[j]),]
			output[k,1:7]<-tmp2[,1:7][1,]
			values<-as.numeric(tmp2$F635.diff.B635)

			output$Intensity1[k]<-values[1]
			output$Intensity2[k]<-values[2]

			k=k+1
		}
	}

	output$Mean<-apply(cbind.data.frame(output$Intensity1,output$Intensity2),1,mean)
	output$SDev<-apply(cbind.data.frame(output$Intensity1,output$Intensity2),1,sd)

	valid.values<-ifelse(output$Intensity1>1 & output$Intensity2>1,TRUE,FALSE)

	output$log1<-ifelse(valid.values, log2(output$Intensity1), 0)
	output$log2<-ifelse(valid.values, log2(output$Intensity2), 0)
	output$log.mean<-ifelse(valid.values, log2(output$Mean), 0)
	output$mean.log<-apply(cbind.data.frame(output$log1,output$log2),1,mean)
	rownames(output)<-seq(1,dim(output)[1])
	
	source("D:\\Stats\\Doc R\\Scripts R\\Kinome\\kinome.loess.R")
	KL<-kinome.loess(output, span=span, graph.title=graph.title)
	
	output$fit<-KL$fit
	output$resid<-KL$resid
	
	output
}