kinome.analysis.v4<-function(X,Y,min.fc=2,min.p=0.05)

{
	
	ref<-cbind.data.frame(log=X$fit,v=X$resid^2,flags=X$Flags)
	exp<-cbind.data.frame(log=Y$fit,v=Y$resid^2,flags=Y$Flags)
	slide<-X$Slide
	wells<-as.character(X$Spot.Coord)
	kinases<-X$ID

	flags<-ifelse(ref$flags=="Flag" | exp$flags=="Flag","Flag","Ok")
	
	# substrates<-X$Substrate
	
	lev<-length(wells)
	init<-rep(0,lev)	

	result<-cbind.data.frame(Slide=slide, well=wells, kinase=kinases,
						log.Ref=round(ref$log,3), log.Exp=round(exp$log,3), Log.R=init, foldchange=init,
						t=init, p.value=init, adj.p=init, status=0, Flags=flags)	#Substrate=substrates,

	Log.R<-round(exp$log-ref$log,3); result$Log.R<-Log.R
	ratio<-2^Log.R; 
	result$foldchange<-ifelse(ratio>=1,ratio,(-1/ratio))
	result$foldchange<-round(result$foldchange,2)
	t<-Log.R/sqrt(ref$v+exp$v)
	result$t<-round(t,3)

	n<-dim(result)[1]
	p.value<-result$p.value

	for (i in 1:n)
		p.value[i]<-(1-pt(abs(t[i]),df=3))*2
	
	result$p.value<-signif(p.value,2)
			
	result$adj.p<-adj.p<-round(p.adjust(p.value,method="BH"),4)
	
	fc<-result$foldchange
	for (i in 1:n)
		{
			if (!is.na(adj.p[i]) & adj.p[i]<0.05 & abs(fc[i])>min.fc) result$status[i]<-ifelse(Log.R[i]<0,"lower","higher")
			else result$status[i]<-"equiv"
		}
	list.lo<-result[which(result$status=="lower"),]
	list.lo<-list.lo[order(list.lo$foldchange),]
	list.up<-result[which(result$status=="higher"),]
	list.up<-list.up[order(list.up$foldchange,decreasing=TRUE),]

	list(result=result, lo=list.lo, up=list.up)

}