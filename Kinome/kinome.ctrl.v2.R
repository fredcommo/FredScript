kinome.ctrl.v2<-function(X,K.type)

{
	internal.controls<-cbind.data.frame(X[1,1:3],n.replic=0,mean=0,sd=0)
	colnames(internal.controls)[2]<-"ctrl.type"	
	{
	if (K.type=="ST2")
		ctrl.list<-c("neg._control")
	else 
		ctrl.list<-c("neg._control","ZZ_control-Y","control")
	}
	
	n<-length(ctrl.list)
	

	for (i in 1:n)
	{
		tmp<-X[which(X$kinase==ctrl.list[i]),]
		tmp.values<-c(tmp$value1,tmp$value2)
		internal.controls [i,1:3]<-tmp[1,1:3]
		internal.controls [i,4]<-length(tmp.values)
		internal.controls [i,5]<-round(mean(tmp.values),0)
		internal.controls [i,6]<-round(sd(tmp.values),0)
	}
	
	internal.controls$ctrl.type<-ctrl.list
	internal.controls
}