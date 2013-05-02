
X<-lame1$F635.diff.B635

neg.deduc<-function(X)
{
	neg.value<-mean(X$F635.diff.B635[which(X$ID=="neg._control")])
	X$F635.diff.B635<-X$F635.diff.B635-neg.value
	X
}

lame1<-neg.deduct(lame1)
lame2<-neg.deduct(lame2)
lame3<-neg.deduct(lame3)
lame4<-neg.deduct(lame4)
lame5<-neg.deduct(lame5)
lame6<-neg.deduct(lame6)


estim<-function(X)
{
	n<-dim(X)[1]
	index<-seq(1,n)

	Slide<-X$Slide[which(index%%2==1)]
	Well<-X$Spot.Coord[which(index%%2==1)]
	Kinase<-X$ID[which(index%%2==1)]

	#data<-X$F635.diff.B635
	
	int.cor1<-X$F635.diff.B635[which(index%%2==1)]
	int.cor2<-X$F635.diff.B635[which(index%%2==0)]

	mu<-rowMeans(cbind(int.cor1,int.cor2))

	x1<-log2(int.cor1)
	x2<-log2(int.cor2)

	x1<-ifelse(int.cor1<0 | int.cor2<0, 0, x1)
	x2<-ifelse(int.cor1<0 | int.cor2<0, 0, x2)

	mu.log<-rowMeans(cbind(x1,x2))
	log.mu<-log2(rowMeans(cbind(int.cor1,int.cor2)))

	X1<-c(x1,x2)
	X2<-c(x2,x1)

	loess.res<-loess(x2~x1)
	new.x<-(x1+x2)/2
	mean.fit<-predict(loess.res,newdata=new.x,se=TRUE)

	na.index<-which(is.na(mean.fit$fit))
		if (length(na.index>0)) 
			{
			mean.fit$fit[na.index]<-loess.res$fitted[na.index]
			mean.fit$se[na.index]<-loess.res$residuals[na.index]
			}

	newdata<-cbind.data.frame(Slide=Slide, Well=Well, Kinase=Kinase,
						int.cor1=int.cor1, int.cor2=int.cor2, mean=mu, log1=x1, log2=x2,
						mean.log=mu.log, log.means=log.mu, fit=mean.fit$fit, resid=mean.fit$se)
	newdata
}

#Res<-cbind.data.frame(x1,x2,fit1=loess.res$fitted[1:408],fit2=loess.res$fitted[409:816],pred=predict(loess.res,newdata=(x1+x2)/2))


# Normalisation

	Slide.type<-levels(Sensib$Slide)

	for (i in 1:length(Slide.type))
	{
	Ctrl.1<-Sensib$F635.Mean[which(Sensib$ID=="control" & Sensib$Slide==Slide.type[i])]
	Ctrl.2<-Resist$F635.Mean[which(Resist$ID=="control" & Resist$Slide==Slide.type[i])]

	cor.fact<-mean(Ctrl.1)/mean(Ctrl.2)
	
	values<-Resist$F635.diff.B635[which(Resist$Slide==Slide.type[i])]
	values<-values*cor.fact
	Resist$F635.diff.B635[which(Resist$Slide==Slide.type[i])]<-values
	}

newS<-estim(Sensib)
newR<-estim(Resist)



diff<-newR$fit-newS$fit
Sdev<-sqrt(newS$resid^2+newR$resid^2)
t<-diff/Sdev
seuil.t<-qt(0.975,1)
p.value<-1-pt(abs(t),df=1)
adj.p<-p.adjust(p.value,method="frd")
ratio<-2^diff
fc<-ifelse(ratio<1, -1/ratio, ratio)

estim1<-2^newS$fit
estim2<-2^newR$fit

Results<-cbind.data.frame(newS[,1:3],estim1,estim2,diff,ratio,fc,t,p.value,adj.p)

filtre<-function(Results,min.p,min.fc,min.int)
{
	value1<-X$estim1
	value2<-X$estim2
	fc<-X$fc
	p.value<-X$p.value

	index.up<-which(fc>=min.fc & p.value<=min.p & (value1>min.int | value2>min.int))
	index.low<-which(fc<=(-min.fc) & p.value<=min.p & (value1>min.int | value2>min.int))

	list(Results, up=Results[index.up,],low=Results[index.low,])
}

color.pt<-function(p.value,fc,min.p=0.05,min.fc=2)
{
	
	n<-length(p.value)
	color<-rep("lightgrey",n)
	for (i in 1:n)
	{
	if (p.value[i]>=min.p) color[i]="lightgrey"
	else
		{
		if(fc[i]>min.fc) color[i]="red"
		if(fc[i]<(-min.fc)) color[i]="green3"
		}
	}
	color
}

plot(-log10(p.value)~log10(ratio),col=color.pt(p.value,fc),pch=19)

MA.plot<-function(X,Y,min.p=0.05,min.fc=2)
{
	M<-Y$fit-X$fit
	A<-log2((X$mean+Y$mean)/2)

	colors<-color.pt(Results$p.value,Results$fc,min.p=min.p,min.fc=min.fc)
	plot(M~A,pch=19,col=colors,xlim=range(0,max(A,na.rm=T)))
}


	R<-cbind.data.frame(int.1,int.2,sd.int.1,sd.int.2,BG.1,BG.2,sd.BG.1,sd.BG.2)

	plot(new2$fit~new3$fit)
	point.cors(new2$fit[index.up]~new3$fit[index.up],pch=19,cex=0.75,col="red")
	point.cors(new2$fit[index.low]~new3$fit[index.low],pch=19,cex=0.75,col="green")

