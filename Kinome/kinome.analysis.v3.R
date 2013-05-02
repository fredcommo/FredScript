kinome.analysis.v3<-function(X,Y)

{
	ref<-cbind.data.frame(log=X$mu.log,v=X$var.log)
	exp<-cbind.data.frame(log=Y$mu.log,v=Y$var.log)
	wells<-X$well
	kinases<-X$kinase
	substrates<-X$substrate
	
	lev<-length(wells)
	init<-rep(0,lev)	

	result<-cbind.data.frame(well=wells,kinase=kinases,
						Substrate=substrates,
						log.Ref=ref$log, log.Exp=exp$log,ratio=init,
						Log.R=init,foldchange=init,t=init,
						p.value=init,adj.p=init,status=0)

	Log.R<-round(exp$log-ref$log,3); result$Log.R<-Log.R
	ratio<-10^Log.R; result$ratio<-ratio
	result$foldchange<-ifelse(ratio>=1,ratio,(-1/ratio))
	t<-Log.R/sqrt(ref$v+exp$v);	result$t<-t

	n<-dim(result)[1]
	p.value<-result$p.value

	for (i in 1:n)
		p.value[i]<-(1-pt(abs(t[i]),df=8))*2
	
	result$p.value<-p.value
			
	result$adj.p<-adj.p<-round(p.adjust(p.value,method="BH"),4)
	
	for (i in 1:n)
		{
			if (!is.na(adj.p[i]) & adj.p[i]<0.05) result$status[i]<-ifelse(Log.R[i]<0,"lower","higher")
			else result$status[i]<-"equiv"
		}
	result

}