FitCurve.nls<-function(x, y, Plot = TRUE, AvoidNeg = TRUE, graph.title=""){

# FitCurve Analysis
# require(MASS)

	x<-as.numeric(x)
	y<-as.numeric(y)

	X.offset<- min(x)

# Smoothing function
	nlm.fit3<-function(bottom, top, k, X.offset, x){
		y.fit<- bottom+(top-bottom)*(1-exp(-k*(x-X.offset)))	#
		return(y.fit)
	}

# formula defined for nls()
	exp.model<-formula(y~bottom+(top-bottom)*(1-exp(-k*(x-X.offset))))

	ks<-ksmooth(x, y, bandwidth=12, kernel="normal")# kernel=15
	X<-ks$x
	Y<-ks$y	
	
	X.lm<-cbind(1,X)
	Y.lm<-Y
	lm.estim<- solve(t(X.lm)%*%X.lm)%*%t(X.lm)%*%Y.lm

	bottom<-Y[1]		# first value
	top<-Y[length(Y)]		# last value
	k=0.1	
		
		
	test<-try(nls(exp.model, data=list(x=X,y=Y), start = list(bottom=bottom, top=top, k=k), algorithm="port"), silent=TRUE)

	{	
	if (class(test)!="try-error"  & lm.estim[2]>0){		# nls is able to fit a exponetial model

		best<-nls(exp.model, data=list(x=X,y=Y), start = list(bottom=bottom, top=top, k=k), algorithm="port")
		bottom=summary(best)$parameters[1]
		top=summary(best)$parameters[2]
		k=summary(best)$parameters[3]

		Y.best<- nlm.fit3(bottom, top, k, X.offset, X)
		choice="Exp.Model"
		Vini<-k*(top-bottom)		
		Vini.fit<-Vini*(X-X.offset) + bottom
		}

	else	{								# the nls-try returned an error
		choice="linear.Model"
		
		rlm.model<-rlm(Y~X, method="MM")

		Y.best<-rlm.model$fitted.values
		a<-rlm.model$coefficients[1]
		b<-rlm.model$coefficients[2]
		Y.best<-a + b*(X)
		Vini.fit<-a + b*(X) 
		Vini<-b
		}
	}

# Convert Neg values to 0 (optional)
	if(AvoidNeg) Vini<-ifelse(Vini<0,1e-3,Vini)

# graphics (optional)
	if (Plot){
		plot(y~x, pch=19, col="blue3", cex=0.5, main=graph.title)
		lines(Y.best~X, col="red")
		lines(Vini.fit~X, col="blue3")
		legend("topleft", legend=choice, bty="n", cex=0.75)
		legend("bottomright", legend=paste("Vini =", round(Vini,3)), bty="n", cex=0.75)
		}	

	return(Vini)

}
# End Analysis
