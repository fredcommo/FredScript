kinome.QC.v3<-function(X,neg.ctrl="neg._control",K.type="ST2",neg.suppr=TRUE,pt.label=TRUE)

# X			: kinases info and duplicates.
# neg.control 	: wells which have to be considered as negative contol (threshold)
# K.type		: Define the type of kinase on the array. 0ptions are "tyr","ST1","ST2" ("ST2" is default)
# neg.suppr		: if TRUE (default), negative values (lower than control) are replaced by NA.
# pt.label		: if TRUE (default) outliers in linear regression are labeled on the plot.

# OUTPUTS
# $X		: values of duplicates, mean, sd ans "flags" (outliers) for the 384 spots (+spots infos)
# $Summary.Duplic : Summary of the linear regression of the duplicates.
# linear regression plot

{
	n<-dim(X)[1]	

	# Appelle kinome.intern.ctrl.R

		source("D:\\Stats\\Doc R\\Scripts R\\Kinome\\kinome.ctrl.v2.R")
		internal.controls<-kinome.ctrl.v2(X,K.type)


	# Appelle kinome.BlankCor.v4

		source("D:\\Stats\\Doc R\\Scripts R\\Kinome\\kinome.BlankCor.v4.R")
		result<-kinome.BlankCor.v4(X,neg.ctrl,neg.suppr)
	

	# Appelle kinome.lm() pour analyse de la régression et sortie graphique

		source("D:\\Stats\\Doc R\\Scripts R\\Kinome\\kinome.lm.v2.R")
		QC<-kinome.lm.v2(result)
		result<-cbind(result,QC=QC)
	
	
	# option identification des points hors IC sur le graphe
		replic1<-result$correct1
		replic2<-result$correct2

		if (pt.label)	
		{
		for (i in 1:n)
			if(!is.na(QC[i]) & QC[i]=="flag")
			{
				points(x=replic1[i], y=replic2[i], pch=20, cex=0.75, col="red")
				legend(x=replic1[i], y=replic2[i], legend=X$well[i], cex=0.5, bty="n",
						x.intersp =-0.5, y.intersp=0.5, xjust=0,yjust=0)
			}
		}

	# output (list : duplicates, regression on duplicates)
		
	
		list(X=result, Internal.Ctrl=internal.controls,Summary.Duplic=summary(lm(replic2~replic1)))

}
