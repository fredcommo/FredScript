FitCurve.nlm6<-function(x, y, tol=1e-12, Plot=TRUE, graph.title=NULL){

	x<-as.numeric(x)
	y<-as.numeric(y)

	X.offset<- min(x)

# Smoothing function
	nlm.fit3<-function(bottom, top, k, X.offset, x){
		y.fit<- bottom+(top-bottom)*(1-exp(-k*(x-X.offset)))
		return(y.fit)
	}

# formula defined for nls()
	exp.model<-formula(y~bottom+(top-bottom)*(1-exp(-k*(x-X.offset))))


# Estimating k (not used)
	k.estim<-function(X, Y){
		X.lm<-cbind(1,c(min(X),max(X)))		
		lm.estim<- solve(t(X.lm)%*%X.lm)%*%t(X.lm)%*%c(Y[1],Y[length(Y)])
		a.temp<-lm.estim[1]
		b.temp<-lm.estim[2]
		y.target<-2/3*(1-1/exp(1))*max(Y)
		x.mid<-(y.target-a.temp)/b.temp
		k<-1/x.mid
		return(k)
	}


# Initialisation
	ks<-ksmooth(x, y, bandwidth=15, n.points=25, kernel="normal")# kernel=15
	X<-ks$x
	Y<-ks$y

	#bottom<-min(Y)
	#top<-max(Y)
	#x.lm<-cbind(1,x[1:10])
	#lm.estim <- solve(t(x.lm)%*%x.lm)%*%t(x.lm)%*%Y[1:10]
	#k = lm.estim[2]/top
	

# model testing
	rlm1<-rlm(y~x, method="MM")
	# rlm2<-rlm(y~I(bottom+(top-bottom)*(1-exp(-k*(x-X.offset)))), method="MM")

	x.lm<-cbind(1,x[1:4])
	lm.estim<- solve(t(x.lm)%*%x.lm)%*%t(x.lm)%*%y[1:4]
	k = lm.estim[2]/top*1.25
	
	{
	if (lm.estim[2]>0.2){

		# initial weights
			rlm.fit<-rlm2$fitted.values
			bottom<-min(rlm.fit)/2
			top<-max(rlm.fit)*2
			weight<-rlm2$w		
			k = lm.estim[2]/top*1.25
		
		# Optimisation
	
			best<-nls(exp.model, data=list(x=x,y=y), start = list(bottom=bottom, top=top, k=k), weights=weight, algorithm="port")
					#, control=nls.control(maxiter=100, tol=tol, minFactor=1e-10))
			
			bottom=summary(best)$parameters[1]
			top=summary(best)$parameters[2]
			k=summary(best)$parameters[3]

			newX<-seq(min(x),max(x),length=100)
			Y.best<- nlm.fit3(bottom, top, k, X.offset, newX)
			Vini.fit<-k*top*(newX-X.offset) + bottom
			Vini<-k*top			
		}

	else{
		a<-rlm1$coefficients[1]
		b<-rlm1$coefficients[2]
		newX<-seq(min(x),max(x),length=100)
		Y.best<-a + b*(newX) 				#(newX - X.offset)
		Vini.fit<-a + b*(newX) 
		Vini<-b				#(newX - X.offset)
		}
	}

# Plot (optional)
	if (Plot){
		plot(y~x, pch=19, col= "blue2", xlab= "Cycle", ylab= "Fluo", main=graph.title)
		lines(Y.best~newX,col= "lightblue3")
		lines(test$fitted.values~x, lty=3)
		lines(Vini.fit~newX, col= "blue3")
		legend("topleft", legend=paste("Vi =", round(Vini,3)), cex=1, bty="n")
	}
	
	return(as.numeric(Vini))
}



