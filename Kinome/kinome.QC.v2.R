kinome.QC.v2<-function(X,pt.label=TRUE)

# X		: corrected values (with Blank deduction).
# pt.label	: if TRUE (default) outliers in linear regression are labeled on the plot.

# OUTPUTS
# $X		: values of duplicates, mean, sd ans "flags" (outliers) for the 384 spots (+spots infos)
# $Summary.Duplic : Summary of the linear regression of the duplicates.
# linear regression plot

{
	X<-X$Values
	n<-dim(X)[1]	
	
	# create a data frame to identify wells
	# data.frame info sequences : fichier de sortie

	result<-data.frame(well=NA,kinase=NA,substrate=NA,
					value1=0,value2=0,mean=0,sdev=0)
			
		i=1	
		for (j in 1:16)
			for (k in 1:24)
			{
				result[i,1]<-paste(LETTERS[j],k,sep="")	# rename wells
				i=i+1
			}	

		n<-dim(result)[1]
		
		# calcul des moyennes et sd par duplicate
		for (i in 1:n)
		{
			well<-result$well[i]
			tmp<-X[which(X$well==well),]
			result$kinase[i]<-as.character(tmp$kinase[1])
			result$substrate[i]<-as.character(tmp$substrate[1])
			result$mean[i]<-mean(tmp$value)
			{
			if (!is.na(result$mean[i]))	# is.na(mean)=TRUE indique qu'au moins 1 des 2 valeurs est NA
				{
				result$value1[i]<-as.numeric(tmp$value[1])
				result$value2[i]<-as.numeric(tmp$value[2])
				result$sdev[i]<-sd(tmp$value)
				}
			else
				{
				result$value1[i]<-NA
				result$value2[i]<-NA
				result$sdev[i]<-NA
				}
			}
		}	

	# Appelle kinome.lm() pour analyse de la régression et sortie graphique

		source("D:\\Stats\\Doc R\\Scripts R\\Kinome\\kinome.lm.R")
		QC<-kinome.lm(result)
		result<-cbind.data.frame(result,QC=QC)

	# option identification des points hors IC sur le graphe
		if (pt.label)	
		{
		for (i in 1:n)
			if(!is.na(result$QC[i]) & result$QC[i]=="flag")
			{
				legend(x=result$value1[i],
						y=result$value2[i],
						legend=result$well[i],
						cex=0.5,bty="n",
						x.intersp =-0.5,y.intersp=0.5,
						xjust=0,yjust=0)
			}
		}

	# output (list : duplicates, regression on duplicates)
	
		replicate1<-result$value1
		replicate2<-result$value2

		list(Values=result, Summary.Duplic=summary(lm(replicate2~replicate1)))

}
