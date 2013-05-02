kinome.BlankCor.v5<-function(X,K.type,neg.ctrl)

	# INPUTS
# X 			: file with wells (column1) and values (column2)
# K.type		: The type of the studied kinase. Options are "tyr", "ST1", "ST2" (default)
# neg.control 	: wells which have to be considered as negative contol (threshold)

	# OUTPUTS
# $Values			: kinases infos, corrected values (duplicates), means, var.
# $intern.ctrl		: values for control wells (number of controls depends on the type of slide)

# begin
	{

		source("D:\\Stats\\Doc R\\Scripts R\\Kinome\\kinome.ctrl.v2.R")
		intern.ctrl<-kinome.ctrl.v2(X,K.type)
	
		# background values and statistics
			bg<-intern.ctrl[which(intern.ctrl$ctrl.type==neg.ctrl),]
			n.bg<-bg$n.replic
			mu.bg<-bg$mean
			var.bg<-bg$sd^2

			t<-qt(0.975,n.bg-1)
			Blank<-mu.bg+t*bg$sd*sqrt(1/n.bg)
	


# outputs :

	results<-cbind.data.frame(X[,1:4], value1=0, value2=0, mean=0,sdev=0)
	
		l<-dim(results)[1]

		for (i in 1:l)
		{
			
			results$value1[i]<-value1<-X$value1[i]-Blank
			results$value2[i]<-value2<-X$value2[i]-Blank
					
			tmp<-c(value1,value2)
			n.i<-length(tmp)
			results$mean[i]<-mu<-round(mean(tmp,na.rm=TRUE),2)
			results$sdev[i]<-S<-round(sqrt(var(tmp,na.rm=TRUE)/n.i+var.bg/n.bg),2)
			
		}

		list(Values=results, intern.ctrl=intern.ctrl)

	} # end