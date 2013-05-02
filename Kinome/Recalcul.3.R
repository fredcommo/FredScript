Recalcul.3<-function(X)

# mu.Log ans var.log calculation after normalization.
# The function also redefine the QC status

{
	new<-cbind.data.frame(X[,1:2],F.cor=0,F.var=0,log.mu=0,log.var=0)

	new$F.cor<-rowMeans(cbind(X[,5],X[,6]))
	new$F.var<-1/4*(X[,7]^2/X[,17]+X[,8]^2/X[,18])
	B.var<-1/4*(X[,13]^2/X[,19]+X[,14]^2/X[,20])

	new$F.var<-new$F.var/2+B.var/2

	new$log.var<-log10(new$F.var/new$F.cor^2+1)
	new$log.mu<-ifelse(new$F.cor<0,1,log10(new$F.cor)-1/2*new$log.var)

	new
}

