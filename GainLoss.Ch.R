GainLoss.Ch<-function(x,...)
{
	index.xname<-which(clones.info.Band$Arm=="xname")
	ep.aCGH.xname<-ep.Band.aCGH[index.xname,keep=T]
	summarize.clones(ep.aCGH.xname)->Summary.xname

	fga.func(ep.aCGH.xname)->fga.xname

	cbind(Sample=as.vector(phenotype(ep.aCGH)$Id),DR=as.vector(phenotype(ep.aCGH)$D.R),gainP=signif(ep.fga$gainP*100,4),lossP=signif(ep.fga$lossP*100,4))->GainLoss.xname

	as.data.frame(GainLoss.xname)->GainLoss.xname
	
	GainLoss.xname
}

