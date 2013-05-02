kinome.BlankCor.v3<-function(X,neg.ctrl="neg._control",K.type="ST2",neg.suppr=FALSE)

	# INPUTS
# X 			: file with wells (column1) and values (column2)
# neg.control 	: wells which have to be considered as negative contol (threshold)
# K.type		: Define the type of kinase on the array. 0ptions are "tyr","ST1","ST2" ("ST2" is default)
# neg.suppr		: if TRUE , negative values (lower than control) are replaced by NA. Default is "FALSE"

	# OUTPUTS
# $Values			: corrected values.
# $internal.Ctrls	: mean and sd values for controls.
 
{


	source("D:\\Stats\\Doc R\\Scripts R\\Kinome\\kinome.ctrl.R")
	internal.controls<-kinome.ctrl(X,K.type=K.type)

		result<-X
		n<-dim(result)[1]
		cor.values<-X[which(X$kinase==neg.ctrl),]$value

		m<-mean(cor.values)
		l<-length(cor.values)
		s<-sd(cor.values)
		t<-qt(0.975,l-1)
		IC<-t*s*sqrt(1/l)
	
		Blank<-m+IC
		
		for (i in 1:n)
			result$value[i]<-result$value[i]-Blank
	
		result$well<-as.factor(result$well)

		if (neg.suppr)
			result$value[which(result$value<=0)]<-NA
	
		list(Values=result,internal.Ctrls=internal.controls)
}

