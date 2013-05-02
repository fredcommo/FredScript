MA.plot<-function(X)
{

	data<-X
	n<-dim(data)[1]
	index<-seq(1,n)

	data$F635<-ifelse(data$F635<=1,1,data$F635)

	val1<-data$F635[which(index%%2==1)]
	val2<-data$F635[which(index%%2==0)]

	x<-log10((val1+val2)/2)
	y<-log10(val1/val2)

	index.up<-which(val1/val2>2)
	index.low<-which(val1/val2<0.5)

	plot(y~x,pch=19,cex=0.75)
	points(y[index.up]~x[index.up],pch=19,cex=0.75,col="red")
	points(y[index.low]~x[index.low],pch=19,cex=0.75,col="green")
	# identify(y~x)
	output<-cbind.data.frame(X[which(index%%2==1),1:6],val1=val1,val2=val2)
	output	
}
