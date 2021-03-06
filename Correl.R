## Calcule et affcihe les corr�lations entre 2 tables

Correl<-function(x,y,..)

{

	data1<-x ;
	data2<-y ;
	
	# concat les donn�es exp�rimentation
	data<-cbind(data1[,-c(1:6)],data2[,-c(1:6)]);
		
	# cr�e un tableau des coef R� adjust.
	coef<-table(colnames(dir[,-c(1:6)]),"coef"=rep(0,17));
	

	#par(mfrow=c(5,4)) ;
	
	for (i in 1:17)
		{
		reg<-lm(data[,i]~data[,i+17],data);
		rownames(coef)[i]<-colnames(data[i]);
		coef[i]<-round(summary(reg)$adj.r.squared,digits=3);
		#plot(data[,i+17]~data[,i],main=colnames(data[i]),xlab="LogR_1",ylab="LogR_2");
		#abline(reg,col="red");
		#legend("bottomright",legend=round(summary(reg)$adj.r.squared,digits=3),bty="n",text.col="red");
		}
	list(coef=coef);
}