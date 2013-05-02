kinome.analysis.v5<-function(X,Y,min.p=0.05,min.fc=2,min.int=25)
{
	diff<-Y$fit-X$fit
	sd1<-as.numeric(X$resid)
	sd2<-as.numeric(Y$resid)
	Sdev<-sqrt(sd1^2+sd2^2)
	t<-diff/Sdev

	seuil.t<-qt(0.975,1)
	p.value<-1-pt(abs(t),df=1)
	adj.p<-p.adjust(p.value,method="BH")

	ratio<-2^diff
	fc<-ifelse(ratio<1, -1/ratio, ratio)

	estim1<-2^X$fit
	estim2<-2^Y$fit

	Results<-cbind.data.frame(X[,1:7],estim1,estim2,diff,ratio,fc,t,p.value,adj.p)

	filtre<-function(Results,min.p,min.fc,min.int)
	{
		value1<-Results$estim1
		value2<-Results$estim2
		fc<-Results$fc
		p.value<-Results$p.value
	
		index.up<-which(fc>=min.fc & p.value<=min.p & (value1>min.int | value2>min.int))
		index.low<-which(fc<=(-min.fc) & p.value<=min.p & (value1>min.int | value2>min.int))
	
		list(Results=Results, up=Results[index.up,],low=Results[index.low,])
	}
	
	output<-filtre(Results, min.p, min.fc, min.int)

	color.pt<-function(Results, index.up, index.low)
	{
	
		n<-dim(Results)[1]
		color<-rep("lightgrey",n)

		index.up<-rownames(index.up)
		index.low<-rownames(index.low)

		for (i in 1:n)
		{
			if(i %in% index.up) color[i]="red"
			if(i %in% index.low) color[i]="green3"
		}
		color
	}


	# MA.plot<-function(Results, index.up, index.low)

		M<-Y$fit-X$fit
		A<-log2((X$Mean+Y$Mean)/2)

		colors<-color.pt(Results, output$up, output$low)
		plot(M~A,pch=19,col=colors,xlim=range(0,max(A,na.rm=T)))

	output

}