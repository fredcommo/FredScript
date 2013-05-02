kinome.test2<-function(X,Y)

# X    	: reference condition (untreated)
# Y	  	: experimental condition

{
	
	ref<-cbind.data.frame(ref.log=X$log.mu,ref.var=X$log.var)
	exp<-cbind.data.frame(exp.log=Y$log.mu,exp.var=Y$log.var)
	wells<-X$well
	kinases<-X$kinase
	#substrates<-X$substrate

		
	lev<-length(wells)
	init<-rep(0,lev)	

	# Fichier de sortie
	result<-cbind.data.frame(well=wells, kinase=kinases,Log.Ref=ref$ref.log, Log.Exp=exp$exp.log, Log.R=init,
						foldchange=init, t=init, p.value=init, adj.p=init)
	
	

	# Remplacement des valeurs négatives (renommées NA) par 1 pour calcul Log
	ref$ref.log[which(is.na(ref$ref.log))]<-1
	exp$exp.log[which(is.na(exp$exp.log))]<-1

	
	diff=0
	ratio=0
	
	
		for (i in 1:lev)
		{
			result$Log.R[i]<-Log.R<-exp$exp.log[i]-ref$ref.log[i]		
			ratio<-10^Log.R

			result$foldchange[i]<-ifelse(ratio<1,round(-1/ratio,0),round(ratio,0))
			result$t[i]<-Log.R/sqrt(exp$exp.var[i]+ref$ref.var[i])
			result$p.value[i]<-pt(abs(result$t[i]),8,lower.tail=FALSE)*2
		}
	
		
		result$adj.p<-round(p.adjust(result$p.value,method="BH"),4)
		
		result
}
