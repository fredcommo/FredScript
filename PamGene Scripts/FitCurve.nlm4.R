FitCurve.nlm4<-function(x, y, tol=1e-12, Plot=TRUE, graph.title=NULL){

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

# Initializing parameters

	lm.model<-lm(y~x)
	a<-lm.model$coefficients[1]
	b<-lm.model$coefficients[2]
	lm.sce<-summary(lm.model)$sigma
	
	Vini<-b


	{
	if (b>0.2){

		# Initialisation	
			bottom<-min(y)
			top<-max(y)
			x.lm<-cbind(1,x)
			lm.estim<- solve(t(x.lm)%*%x.lm)%*%t(x.lm)%*%y
			k = lm.estim[2]/top
			
		# initial weights
			test<-rlm(y~I(bottom+(top-bottom)*(1-exp(-k*(x-X.offset)))))
			weight<-test$w		
			ks<-ksmooth(x, test$fitted.values, bandwidth=12, kernel="normal")		# kernel=15
			X<-ks$x
			Y<-ks$y

			
			

		# Optimisation
			best<-nls(exp.model, data=list(x=c(x,x),y=c(y,test$fitted.values)), start = list(bottom=bottom, top=top, k=k), weights=weight,
					algorithm="port", control=nls.control(tol=tol))
			bottom=summary(best)$parameters[1]
			top=summary(best)$parameters[2]
			k=summary(best)$parameters[3]

			Y.fit<- nlm.fit3(bottom, top, k, X.offset, x)
			exp.sce<-sqrt(sum(y-Y.fit)^2/(length(x)-3))
	

		# Fitting values: the choise between lm and exp.model depends on best residuals
			{
			if(exp.sce<lm.sce){
				Vini<-k*top
				newX<-seq(min(x), max(x), length=100)
				Y.best<- nlm.fit3(bottom, top, k, X.offset, newX)
				Vini.fit<-k*top*(newX-X.offset) + bottom	#+Y.best[1]
				}
			else{
				newX<-seq(min(x), max(x), length=100)
				Y.best<-a + b*(newX) 			#(newX - X.offset)
				Vini.fit<-a + b*(newX) 			#(newX - X.offset)
				}
			}
		}

	else{
		newX<-seq(min(x),max(x),length=100)
		Y.best<-a + b*(newX) 				#(newX - X.offset)
		Vini.fit<-a + b*(newX) 				#(newX - X.offset)
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



