FitCurve.nlm<-function(x, y, tol=1e-12, Plot=TRUE, graph.title=NULL){

	x<-as.numeric(x)
	y<-as.numeric(y)

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


	X<-c(x[1], x[1]+(x[2]-x[1])*0.3, x[1]+(x[2]-x[1])*0.6, x[-1])
	Y<-c(y[1], y[1]+(y[2]-y[1])*0.35, y[1]+(y[2]-y[1])*0.70, y[-1])


	#plot(Y~X, pch=19, cex=1, col="blue2")
	#k=0.1

	x.lm<-cbind(1,x)
	lm.estim<- solve(t(x.lm)%*%x.lm)%*%t(x.lm)%*%y
	
	{
	if (lm.estim[2]>0.2){

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
		Vini<-k*top

	# Fitting values
		newX<-seq(min(X),max(X),length=100)
		Y.best<- nlm.fit3(bottom, top, k, X.offset, newX)
		Vini.fit<-k*top*(X-X.offset)+Y.best[1]
		}
	else{
		a<-Vini<-lm.estim[2]
		b<-lm.estim[1]
		newX<-seq(min(X),max(X),length=100)
		Y.best<-a*(newX-X.offset)+b
		Vini.fit<-a*(X-X.offset)+b
		}
	}

# Plot
	if (Plot){
		plot(y~x, pch=19, col= "blue2", xlab= "Cycle", ylab= "Fluo", main=graph.title)
		lines(Y.best~newX,col= "lightblue3")
		lines(Vini.fit~X, col= "blue3")
		legend("topleft", legend=paste("Vi =", round(Vini,3)), cex=1.25, bty="n")
	}
	
	return(Vini)
}



