FitCurve.nlm2<-function(x, y, tol=1e-12, Plot=TRUE, graph.title=NULL){

	x<-as.numeric(x)
	y<-as.numeric(y)

	X.offset<-min(x)

# Smoothing function
	nlm.fit3<-function(bottom, top, k, X.offset, x){
		y.fit<- bottom+(top-bottom)*(1-exp(-k*(x-X.offset)))
		return(y.fit)
	}


# Sum of Squares
	sce <- function(param, x, yobs) {
		b<-param[1]
		t<-param[2]
		k<-param[3]
		ytheo <- nlm.fit3(bottom=b, top=t, k=k, X.offset, x)
		return(sum((yobs - ytheo)^2))
	}

# Initializing parameters
	#X<-c(x[1]*0.5, x, x[length(x)]*1.25)		pour HillFit()
	#Y<-c(floor(y[1]), y, ceiling(y[length(y)]))	pour HillFit()


	lm.model<-lm(y~x)
	a<-lm.model$coefficients[1]
	b<-lm.model$coefficients[2]
	lm.sce<-summary(lm.model)$sigma

	{
	if (b>0.1){
		
		y2<-y[1]+(y[2]-y[1])*0.35
		y3<-y[1]+(y[2]-y[1])*0.70
		X<-c(x[1], x[1]+(x[2]-x[1])*0.3, x[1]+(x[2]-x[1])*0.6, x[-1])
		Y<-c(y[1], y2, y3, y[-1])

		for (iter in 1:5){
		# Initialisation	
			X.offset<-min(X)
			bottom<-min(Y)
			top<-max(Y)
			X.lm<-cbind(1,X[1:4])
			k.estim<- solve(t(X.lm)%*%X.lm)%*%t(X.lm)%*%Y[1:4]
			k= k.estim[2]/top*1.1

		# Optimisation
			best<-nlm(f = sce, p = c(bottom, top, k), x = X, yobs = Y, steptol=tol)

			bottom<-best$estimate[1]
			top<-best$estimate[2]
			k<-best$estimate[3]
			y2<-nlm.fit3(bottom, top, k, X.offset, X[2])
			y3<-nlm.fit3(bottom, top, k, X.offset, X[3])
			}

		Y.fit<- nlm.fit3(bottom, top, k, X.offset, x)
		exp.sce<-sqrt(sum(y-Y.fit)^2/(length(x)-3))

		# Fitting values
		{
		if(exp.sce<lm.sce){
			Vini<-k*top
			newX<-seq(min(X),max(X),length=100)
			Y.best<- nlm.fit3(bottom, top, k, X.offset, newX)
			Vini.fit<-k*top*(newX-X.offset)+Y.best[1]
			}
		else{
			newX<-seq(min(X),max(X),length=100)
			Y.best<-a+b*(newX-X.offset)
			Vini.fit<-a+b*(newX-X.offset)
			Vini<-b
			}
		}
	}

	else{
		newX<-seq(min(x),max(x),length=100)
		Y.best<-a+b*(newX-X.offset)
		Vini.fit<-a+b*(newX-X.offset)
		Vini<-b
		}
	}

# Plot
	if (Plot){
		plot(y~x, pch=19, col= "blue2", xlab= "Cycle", ylab= "Fluo", main=graph.title)
		lines(Y.best~newX,col= "lightblue3")
		lines(Vini.fit~newX, col= "blue3")
		legend("topleft", legend=paste("Vi =", round(Vini,3)), cex=1, bty="n")
	}
	
	return(as.numeric(Vini))
}



