kinome.ctrl<-function(X,K.type)

{
	internal.controls<-cbind.data.frame(control.type=NA,X[1,c(4,6)],n.replic=0,mean=0,s.dev=0)	
	{
	if (K.type=="ST2")
		ctrl.list<-c("neg._control","Blank")
	else 
		ctrl.list<-c("neg._control","ZZ_control-Y","control","Blank")
	}
	
	n<-length(ctrl.list)
	

	for (i in 1:(n-1))
	{
		tmp<-X[which(X$kinase==ctrl.list[i]),]
		internal.controls [i,2:3]<-tmp[1,c(4,6)]
		internal.controls [i,4]<-length(tmp$value)
		internal.controls [i,5]<-round(mean(tmp$value),0)
		internal.controls [i,6]<-round(sd(tmp$value),0)
	}
	
	tmp<-X[which(X$well=="Blank"),]
		internal.controls [n,2:3]<-tmp[1,c(4,6)]
		internal.controls [n,4]<-length(tmp$value)
		internal.controls [n,5]<-round(mean(tmp$value),0)
		internal.controls [n,6]<-round(sd(tmp$value),)	

	internal.controls$control.type<-ctrl.list
	internal.controls
}