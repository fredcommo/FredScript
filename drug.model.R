drug.model<-function(X,Y,observed,Log=FALSE,title=NULL)

# calcule la surface de réponses en % à partir des coefficients d'un glm

# X 		: drogue1
# Y 		: drogue2
# observed 	: réponses observées
# Log		: transforme les doses drogue en log10(dose). Faux par défaut.

{ 

	obs.logit<-log10(observed/(1-observed))
	
	{
	if (Log)
		{
		drug1<-log10(X)
		drug2<-log10(Y)
		}
	else
		{
		drug1<-X
		drug2<-Y
		}
	}

	model<-glm(obs.logit~drug1*drug2)
	coef<-model$coefficients

	f.model<-function(newX,newY,coef){
	res<-matrix(0,length(newX),length(newY))
	for (i in 1:length(newX))
		{for (j in 1:length(newY))
			{drug<-c(1,newX[i],newY[j],newX[i]*newY[j])
			 res[i,j]<-as.vector(coef)%*%as.vector(drug)
			}
		}
	return(res)
	}

	n1<-length(unique(X))
	n2<-length(unique(Y))
	newX<-seq(0,max(X),len=(n1*4))
	newY<-seq(0,max(Y),len=(n2*4))
	
	if (Log)
		{
			newX<-log10(newX[-1])
			newY<-log10(newY[-1])
		}

	fit.log<-f.model(newX,newY,coef)
	fit.resp<-1/(1+10^(-fit.log))

	persp(fit.resp,theta=35,phi=25,col="lightblue1",
			xlab="drug1",ylab="drug2",zlab="Estim.Resp",main=title)
	
	res.table<-data.frame(coef)
	text=signif(coef[1],3)

	for (i in 2:dim(res.table)[1])
		text<-paste(text ," + ", signif(res.table[i,1],3),".",rownames(res.table)[i])		

	legend(x=-0.60,y=0.39,legend=paste("logit = ",text),bty="n",text.col="blue",cex=1)	

	list(model=summary(model),fit.log=fit.log,fit.resp=fit.resp)
}

