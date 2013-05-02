kinome.BlankCor.v4<-function(X,neg.ctrl,neg.suppr)

	# INPUTS
# X 			: file with wells (column1) and values (column2)
# neg.control 	: wells which have to be considered as negative contol (threshold)
# neg.suppr		: if TRUE (default), negative values (lower than control) are replaced by NA.

	# OUTPUTS
# $X			: means, log.means, var, Log.var.


# begin
	{

		# background values and statistics
			bg1<-X[which(X$kinase==neg.ctrl),]$value1
			bg2<-X[which(X$kinase==neg.ctrl),]$value2

			n.bg<-length(c(bg1,bg2))
			mu.bg<-mean(c(bg1,bg2))
			var.bg<-var(c(bg1,bg2))

# outputs :

	results<-cbind.data.frame(X[,1:4],correct1=0,correct2=0, mean=0,var=0,mu.log=0,var.log=0)
	
		l<-dim(results)[1]

		for (i in 1:l)
		{
			
			correct1<-X$value1[i]-mu.bg
			correct2<-X$value2[i]-mu.bg
		
			
			if (neg.suppr)
			{
				if (correct1<0 | correct2<0) correct1=correct2=NA
			}
			
			tmp<-c(correct1,correct2)
		
			n.i<-length(which(!is.na(tmp)))
			mu<-mean(tmp,na.rm=TRUE)
			
			#{
			#	if(is.na(var(tmp,na.rm=TRUE)))	s2<- var.bg/n.bg
			#	else 	s2<-var(tmp,na.rm=TRUE)/n.i+ var.bg/n.bg
			#}

			s2<-var(tmp,na.rm=TRUE)/n.i+ var.bg/n.bg

			results$correct1[i]<-correct1
			results$correct2[i]<-correct2
			results$mean[i]<-mu
			results$var[i]<-s2
				
			var.log<-log10(s2/mu^2+1)
			mu.log<-log10(mu)-1/2*var.log

			results$mu.log[i]<-mu.log
			results$var.log[i]<-var.log
		}
		results$mu.log[which(is.na(results$mu.log))]<-1
		results$var.log[which(is.na(results$var.log))]<-1
		results

	} # end