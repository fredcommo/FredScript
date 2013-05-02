# FitCurve Analysis
# require(MASS)

# Same as FitCurve.nls, except that there is no kernel pre-smoothing (in PamGene Scripts)

FitCurve.nls2<-function(x, y, Plot = TRUE, AvoidNeg = TRUE, graph.title=""){

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

# Initialization
	bottom<-min(y)		# first value
	top<-max(y)		# last value
	k=0.1	
	
	test<-try(nls(exp.model, data=list(x=x,y=y), start = list(bottom=bottom, top=top, k=k), algorithm="port"), silent=TRUE)

	{	
	if (class(test)!="try-error"){		# is nls able to fit a exponetial model?

		best<-nls(exp.model, data=list(x=x,y=y), start = list(bottom=bottom, top=top, k=k), algorithm="port")
		bottom=summary(best)$parameters[1]
		top=summary(best)$parameters[2]
		k=summary(best)$parameters[3]

		newx<-seq(min(x), max(x), length=100)
		Y.best<- nlm.fit3(bottom, top, k, X.offset, newx)
		choice="Exp.Model"
		Vini<-k*(top-bottom)		
		Vini.fit<-Vini*(x-X.offset) + bottom
		}

	else	{								# the nls-try returned an error
		choice="linear.Model"
		
		rlm.model<-rlm(y~x, method="MM")

		Y.best<-rlm.model$fitted.values
		a<-rlm.model$coefficients[1]
		b<-rlm.model$coefficients[2]
		Y.best<-a + b*(x)
		Vini.fit<-a + b*(x) 
		Vini<-b
		}
	}

# Convert Neg values to 0 (optional)
	if(AvoidNeg) Vini<-ifelse(Vini<0,1e-3,Vini)

# graphics (optional)
	if (Plot){
		plot(y~x, pch=19, col="blue3", cex=0.5, main=graph.title)
		lines(Y.best~newx, col="red")
		lines(Vini.fit~x, col="blue3")
		model = choice
		leg1 <- paste("bottom = ",round(bottom,3), sep="")
		leg2 <- paste("top = ",round(top,3), sep="")
		leg3 <- paste("Vinit = ",round(Vini,3), sep="")
		legend("bottomright", legend=c(model, leg1, leg2, leg3), bty="n", cex=0.75)
		}	

	return(Vini)

}
# End Analysis
