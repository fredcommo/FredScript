HillFit2 <- function (X, Y, test.rlm=TRUE, Slope=TRUE, xlab="", ylab="", title=NULL, print=TRUE){

	# Calcul de l'alignement
		nlm.fit<-function(bottom,top,slope,X50,X){
			Y.fit<-bottom+(top-bottom)/(1+10^((X50-X)*slope))
			return(Y.fit)
			}

	# Fonction sce (somme carré résidus)
		sce <- function(param, X, yobs) {
			bottom<-param[1]
			top<-param[2]
			slope<-param[3]
			X50<-param[4]
			ytheo <- nlm.fit(bottom, top, slope, X50, X)
			return(sum((yobs - ytheo)^2))
			}

	# parametres initiaux
	
		lin.mod <- rlm(Y ~ X)
		init.bottom = min(Y, na.rm = T)
		init.top = Y[order(X)][length(Y)]
		init.slope = coef(lin.mod)[2]
		init.X50 = (max(X, na.rm = T)+min(X, na.rm = T))/2

	plot(Y~X, pch=19, col="blue", main="4P Logistic Regression", sub=title, xlab=xlab, ylab=ylab)
	
	# Y.init<-nlm.fit(init.bottom, init.top, init.slope, init.X50, X)

	cond = TRUE
	if (test.rlm) cond = abs(coef(lin.mod)[2])>0.005

	{
	if (cond) {

		init=c(bottom=init.bottom, top=init.top, slope=init.slope, X50=init.X50)	# initialisation des paramètres
		best<-nlm(f = sce, p = init, X = X, yobs = Y)						# calcul des paramètres
	
		# Récupération des paramètres
			best.bottom <- best$estimate[1]
			best.top <- best$estimate[2]
			Result <- best.slope<-best$estimate[3]
			best.X50 <- best$estimate[4]

		SStot <- length(X)*var(X)
		SSres <- sce(c(best.bottom, best.top, best.slope, best.X50), X, Y)
		r2 = (1-SSres/SStot)

		newX <- seq(min(X), max(X), length=100)						
		Y.best <- nlm.fit(best.bottom, best.top, best.slope, best.X50, newX)
		lines(Y.best~newX, col="orangered", lwd=2)
		
		if (Slope){
			Y50 <- nlm.fit(best.bottom, best.top, best.slope, best.X50, best.X50)
			B <- 10^best.slope
			A <- Y50 - B*best.X50	
			lines ((A + B*newX) ~ newX, lty=3, col="violetred4")
		}

		legend.pos = "topleft"
		if(best.slope<0) legend.pos = "bottomleft"

		legend(legend.pos, legend=c("observed","4PL"), pch=c(19,-1), lty=c(-1,1), lwd=3,
		col=c("blue","orangered"), bty="n")

		}

	else {
		best <- lin.mod
		Result <- coef(lin.mod)[2]
		abline(lin.mod, col="violetred4")
		}
	}
	
	if (print)	print(best)
	Result
}

# 