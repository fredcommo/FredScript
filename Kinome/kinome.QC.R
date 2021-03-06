kinome.QC<-function(X,pt.label=TRUE)
{
	n<-dim(X)[1]	
	
	# QC sur les contr�les
		neg<- X[which(X$kinase=="neg_control"),]
		zz<-X[which(X$kinase=="ZZ_control-Y"),]
		ctrl<- X[which(X$kinase=="control"),]
		blank<-X[which(X$kinase=="Blank"),]

		internal.controls<-data.frame(X[1,c(3,4,6)],n.replic=0,mean=0,s.dev=0)

		internal.controls [1,1:3]<-neg[1,c(3,4,6)]
		internal.controls [1,4]<-length(neg$value)
		internal.controls [1,5]<-mean(neg$value)
		internal.controls [1,6]<-sd(neg$value)
	
		internal.controls [2,1:3]<-zz[1,c(3,4,6)]
		internal.controls [2,4]<-length(zz$value)
		internal.controls [2,5]<-mean(zz$value)
		internal.controls [2,6]<-sd(zz$value)

		internal.controls [3,1:3]<-ctrl[1,c(3,4,6)]
		internal.controls [3,4]<-length(ctrl$value)
		internal.controls [3,5]<-mean(ctrl$value)
		internal.controls [3,6]<-sd(ctrl$value)

		internal.controls [4,1:3]<-blank[1,c(3,4,6)]
		internal.controls [4,4]<-length(blank$value)
		internal.controls [4,5]<-mean(blank$value)
		internal.controls [4,6]<-sd(blank$value)


	# create a data frame to identify wells
	# data.frame info sequences : fichier de sortie

		result<-data.frame(well=X[1,1],kinase=X[1,3],substrate=X[1,5],
					value1=0,value2=0,mean=0,sdev=0)			
		i=1	
		for (j in 1:16)
			for (k in 1:24)
			{
				result[i,1]<-paste(LETTERS[j],k,sep="")	# renomme les wells
				i=i+1
			}	

		n<-dim(result)[1]
		
		# calcul des moyennes et sd par duplicate
		for (i in 1:n)
		{
			well<-result$well[i]
			tmp<-X[which(X$well==well),]
			result$kinase[i]<-tmp$kinase[1]
			result$substrate[i]<-tmp$substrate[1]
			result$value1[i]<-tmp$value[1]
			result$value2[i]<-tmp$value[2]
			result$mean[i]<-mean(tmp$value)
			result$sdev[i]<-sd(tmp$value)
		}	

	# Appelle kinome.lm pour analyse de la r�gression et sortie graphique

		source("D:\\Stats\\Doc R\\Scripts R\\kinome.lm.R")
		QC<-kinome.lm(result)
		result<-cbind.data.frame(result,QC=QC)

	# option identification des points hors IC sur le graphe
		if (pt.label)	
		{
		for (i in 1:n)
			if(result$QC[i]=="flag")
			{
				legend(x=result$value1[i],
						y=result$value2[i],
						legend=result$well[i],
						cex=0.5,bty="n",
						x.intersp =-0.5,y.intersp=0.5,
						xjust=0,yjust=0)
			}
		}

	# sortie valeurs (liste : duplicats, contr�les internes, regression duplicates
	
		replicate1<-result$value1
		replicate2<-result$value2

		list(X=result, internal.Ctrl=internal.controls,
				Summary.Duplic=summary(lm(replicate2~replicate1)))

}
