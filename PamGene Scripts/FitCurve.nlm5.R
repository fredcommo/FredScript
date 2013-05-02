FitCurve.nlm5<-function(x, y, tol=1e-12, Plot=TRUE, graph.title=NULL){

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
	ks<-ksmooth(x, y, bandwidth=12, n.points=25, kernel="normal")# kernel=15
	X<-ks$x
	Y<-ks$y

	bottom<-min(Y)
	top<-max(Y)
	x.lm<-cbind(1,X[1:10])
	lm.estim<- solve(t(x.lm)%*%x.lm)%*%t(x.lm)%*%Y[1:10]
	k = lm.estim[2]/top
	

# model testing
			
	rlm1<-rlm(Y~X, method="MM")
	rlm2<-rlm(Y~I(bottom+(top-bottom)*(1-exp(-k*(X-X.offset)))), method="MM")
	comp<-anova(rlm1, rlm2)
	RSS.ratio<-max(comp$RSS)/min(comp$RSS)
	df1<-length(x)-2
	df2<-length(x)-3
	p<-pf(RSS.ratio, df1,df2, lower.tail=FALSE)
	
	{
	if (p<=0.25 & lm.estim[2]>0){

		# initial weights
		rlm.fit<-rlm2$fitted.values
		bottom<-min(rlm.fit)
		top<-max(rlm.fit)
		weight<-rlm2$w		
		k = lm.estim[2]/top*1.25

		# Optimisation
	
			best<-nls(exp.model, data=list(x=X,y=Y), start = list(bottom=bottom, top=top, k=k), 
					weights=weight, algorithm="port")

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
		lines(rlm1$fitted.values~X, lty=3)
		lines(rlm2$fitted.values~X, lty=3, col="red")
		lines(Vini.fit~newX, col= "blue3")
		legend("topleft", legend=paste("Vi =", round(Vini,3)), cex=1, bty="n")
	}
	
	return(as.numeric(Vini))
}



