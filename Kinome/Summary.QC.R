Summary.QC<-function(X)

# X : QC file from kinome.QC()

{
	x<-X$Values
	n<-dim(x)[1]
	
		pts.OK<-round(length(which(x$QC=="OK"))/n*100,2)
		pts.flag<-round(length(which(x$QC=="flag"))/n*100,2)
		pts.suppr<-round(length(which(x$QC=="suppr"))/n*100,2)

	QC.result<-rbind(pts.OK=pts.OK,pts.flag=pts.flag,pts.suppr=pts.suppr)
	QC.result<-as.data.frame(QC.result)
	colnames(QC.result)<-"Prop"

	list(QC.duplicates=QC.result,Intern.Ctrl=X$Intern.Ctrl,Duplic.reg=X$Summary.Duplic)
}