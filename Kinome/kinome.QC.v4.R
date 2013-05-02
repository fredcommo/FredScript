kinome.QC.v4<-function(X,neg.ctrl="neg._control",K.type="ST2",pt.label=TRUE)

# X			: kinases info and duplicates.
# neg.control 	: wells which have to be considered as negative contol (threshold)
# K.type		: Define the type of kinase on the array. 0ptions are "tyr","ST1","ST2" ("ST2" is default)
# pt.label		: if TRUE (default) outliers in linear regression are labeled on the plot.

# OUTPUTS
# $X		: values of duplicates, mean, sd ans "flags" (outliers) for the 384 spots (+spots infos)
# $Summary.Duplic : Summary of the linear regression of the duplicates.
# linear regression plot

{
	n<-dim(X)[1]	

	
	# Appelle kinome.BlankCor.v5

		source("D:\\Stats\\Doc R\\Scripts R\\Kinome\\kinome.BlankCor.v5.R")
		X.bg<-kinome.BlankCor.v5(X,K.type,neg.ctrl)
	

	# Appelle kinome.lm.v2() pour analyse de la régression et sortie graphique

		#source("D:\\Stats\\Doc R\\Scripts R\\Kinome\\kinome.lm.v2.R")
		#QC<-kinome.lm.v2(X.bg$Values)

		source("D:\\Stats\\Doc R\\Scripts R\\Kinome\\kinome.lm.v3.R")
		QC<-kinome.lm.v3(X.bg$Values)
		X.bg$Values<-cbind(X.bg$Values,QC=QC)

		replic1<-X.bg$Values$value1
		replic2<-X.bg$Values$value2		

	# Marquage des points < BG
		for (i in 1:n)
			if(!is.na(QC[i]) & (QC[i]=="suppr"))
			{
				points(x=replic1[i], y=replic2[i], pch=20, cex=0.75, col="black")
			}
	

	# option identification des points hors IC sur le graphe
		if (pt.label)	
		{
		for (i in 1:n)
			if(!is.na(QC[i]) & (QC[i]=="flag"))
			{
				points(x=replic1[i], y=replic2[i], pch=20, cex=0.75, col="red")
				legend(x=replic1[i], y=replic2[i], legend=X.bg$Values$well[i], cex=0.5, bty="n",
						x.intersp =-0.5, y.intersp=0.5, xjust=0,yjust=0)
			}
		}
	# Marquage de la zone de définition Background
		abline(h=0,lty=3,col="lightblue3")
		abline(v=0,lty=3,col="lightblue3")
		legend(x=min(replic1),y=max(replic2)*3/4,legend="bground",bty="n",cex=0.75,text.col="lightblue4")
		legend("bottomright",legend="bground",bty="n",cex=0.75,text.col="lightblue4")
		
		
	# output (list : BG corrected values, internal ctrls, regression summary on duplicates)
	
		list(Values=X.bg$Values, Intern.Ctrl=X.bg$intern.ctrl,Summary.Duplic=summary(lm(replic2~replic1)))

}
