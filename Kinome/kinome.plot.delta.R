kinome.plot.delta<-function(X, min.p=0.05, min.fc=3, min.int=50, exclude.flags=TRUE)

{
	datas<-X
	
	slide.type<-levels(datas$Slide)

	source("D:\\Stats\\Doc R\\Scripts R\\kinome\\kinome.SlideGrid.R")
	slide.loc<-kinome.SlideGrid()
	
	n<-dim(slide.loc)[1]
	p<-dim(slide.loc)[2]

# Calcul des Scores et représentation

	slide.value<-matrix(0,n,p)
	
	par(mfrow=c(1,4))

	for (k in 1:length(slide.type))
	{
		sub<-datas[which(datas$Slide==slide.type[k]),]

		{
			if (exclude.flags) flags<-sub$Flags
			else flags<-rep("OK",dim(datas)[1])
		}
	
		plot(slide.value,xlim=range(1,p),ylim=range(1:n),type="n",
			xlab="columns",ylab="rows",main=paste("Changes on",slide.type[k]))	

		for (i in 1:n)
			for (j in 1:p)
			{
				well<-slide.loc[i,j]
				
				if(well=="x") points(x=j,y=(54-i),pch=19,cex=0.5,col="lightblue")

				else 
					{
						pval<-sub$p.value[which(sub$Spot.Coord==well)]
						fc<-sub$fc[which(sub$Spot.Coord==well)]
						flag.spot<-flags[which(sub$Spot.Coord==well)]
						int1<-sub$estim1[which(sub$Spot.Coord==well)]
						int2<-sub$estim2[which(sub$Spot.Coord==well)]
						int<-ifelse(int1<min.int & int2<min.int,FALSE,TRUE)

						{
						if (pval>min.p | abs(fc)<min.fc | flag.spot=="Flag" | !int)
							{
								point.color="white"
							}
						
						else 
							#if(pval<=min.p & abs(fc)>min.fc & flag.spot!="Flag")
							{
								if (fc<0) point.color="green"
								else      point.color="red"
							}
	
							points(x=j,y=(54-i),pch=19,cex=log10(abs(fc)),col=point.color)
						}
					}

			}
	}	
}
