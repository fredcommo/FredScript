FitCurve.nls4<-function(x, y, Plot = TRUE, AvoidNeg = TRUE, graph.title=""){

# FitCurve Analysis
# require(MASS)
# use a presmoothing method (kernel)

	X<-as.numeric(x)
	Y<-as.numeric(y)

	X.offset<- min(X, na.rm=T)

# Smoothing function
	nlm.fit3<-function(bottom, top, k, X.offset, x){
		y.fit<- bottom+(top-bottom)*(1-exp(-k*(x-X.offset)))	#
		return(y.fit)
	}

# formula defined for nls()
	exp.model<-formula(y~bottom+(top-bottom)*(1-exp(-k*(x-X.offset))))

			# ks<-ksmooth(x, y, bandwidth=0.25, kernel="normal")	# kernel can be used in the number of values is too small.
	#X <- x	# X<-ks$x
	#Y <- y	# Y<-ks$y	

# Suppress NA values, Inf values, if log has been previously used for Y computation
	
	suppr <- c(which(is.na(Y)), which(Y==(-Inf) | Y==(Inf)))
	if (length(suppr)>0){
		Y <- Y[-suppr]
		X <- X[-suppr]
	}

	new.x <- seq(min(X), max(X), length=50)	# new x.values for graphics
	
	{
	if (length(unique(X))>3){

		# initial evaluation: slope must be >0 to test the exp.model, rlm will be used, otherwise.	
			X.lm<-cbind(1,X)
			Y.lm<-Y
			lm.estim<- solve(t(X.lm)%*%X.lm)%*%t(X.lm)%*%Y.lm

		# initialize the model parameters
			bottom<-max(Y, na.rm=T)					# first value
			top<-min(Y, na.rm=T)					# last value
			k=lm.estim[2]	

	
			test<-try(nls(exp.model, data=list(x=X,y=Y), start = list(bottom=bottom, top=top, k=k)), silent=TRUE)	# nls(...,algorithm="port")
	
			{	
			if (class(test)!="try-error"  & abs(lm.estim[2])>0.001){		# nls is able to fit a exponetial model
	
				best<-nls(exp.model, data=list(x=X,y=Y), start = list(bottom=bottom, top=top, k=k))
				bottom = coef(best)[1]
				top = coef(best)[2]
				k = coef(best)[3]
			
				Y.best<- nlm.fit3(bottom, top, k, X.offset, new.x)
				choice="Exp.Model"
				Vini<-k*(top-bottom)		
				Vini.fit<-Vini*(X-X.offset) + bottom
				}
	
			else	{								# if the nls-try has returned an error
				choice="linear.Model"
			
				rlm.model<-rlm(Y~X, method="M")			# "MM"
	
				Y.best<-rlm.model$fitted.values
				a<-rlm.model$coefficients[1]
				b<-rlm.model$coefficients[2]
				Y.best<-a + b*(new.x)
				Vini.fit<-a + b*(X) 
				Vini<-b
				}
			}
		}
	else	{
		choice="Null slope"
		Y.best=rep(0, length(new.x))
		Vini.fit=rep(0, length(x))
		Vini<-1e-3
		}
	}

# Convert Neg values to 0 (optional)
	if(AvoidNeg) Vini<-ifelse(Vini<0,1e-3,Vini)

# graphics (optional)
	if (Plot & length(unique(X))>3){
		plot(Y~X, pch=19, col="blue3", cex=0.5, main=graph.title)
		lines(Y.best~new.x, col="red")
		lines(Vini.fit~X, col="blue3")
		legend("topleft", legend=choice, bty="n", cex=0.75)
		legend("bottomright", legend=paste("Vini =", round(Vini,3)), bty="n", cex=0.75)
		}	

	return(Vini)

}
# End Analysis
