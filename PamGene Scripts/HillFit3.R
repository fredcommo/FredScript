HillFit3 <- function (X, Y, title=NULL, print=TRUE){

	# Calcul de l'alignement
		nlm.fit<-function(bottom, top, slope, Assym, X50, X){
			Y.fit<-bottom+(top-bottom)/((1+10^((X50-X)*slope))^Assym)
			return(Y.fit)
			}

	# Fonction sce (somme carré résidus)
		sce <- function(param, X, yobs) {
			bottom <- param[1]
			top <- param[2]
			slope <- param[3]
			Assym <- param[4]
			X50 <- param[5]
			ytheo <- nlm.fit(bottom, top, slope, Assym, X50, X)
			return(sum((yobs - ytheo)^2))
			}

	# parametres initiaux
	
		lin.mod <- rlm(Y ~ X)
		init.bottom = min(Y, na.rm = T)
		init.top = max(Y, na.rm = T)
		init.slope = coef(lin.mod)[2]
		init.Assym = 1
		init.X50 = (max(X, na.rm = T)+min(X, na.rm = T))/2

	plot(Y~X, pch=19, col="blue", main="4P Logistic Regression", xlab="X", ylab="Resp (in prop of max)")
	
	# Y.init<-nlm.fit(init.bottom, init.top, init.slope, init.X50, X)


	{
	if (abs(coef(lin.mod)[2])>0.005) {

		init=c(bottom=init.bottom, top=init.top, slope=init.slope, Assym = init.Assym, X50=init.X50)	# initialisation des paramètres
		best<-nlm(f = sce, p = init, X = X, yobs = Y)						# calcul des paramètres
	
		# Récupération des paramètres
			best.bottom <- best$estimate[1]
			best.top <- best$estimate[2]
			Result <- best.slope<-best$estimate[3]
			best.Assym <- best$estimate[4]
			best.X50 <- best$estimate[5]

		SStot <- length(X)*var(X)
		SSres <- sce(c(best.bottom, best.top, best.slope, best.Assym, best.X50), X, Y)
		r2 = (1-SSres/SStot)

		newX <- seq(min(X), max(X), length=100)						
		Y.best <- nlm.fit(best.bottom, best.top, best.slope, best.Assym, best.X50, newX)

		# E50 <- nlm.fit(best.bottom, best.top, best.slope, best.Assym, best.X50, best.X50)

		lines(Y.best~newX, col="orangered", lwd=2)
		#abline (h=E50, v=best.X50, lty=3)
		legend("topleft", legend=c("observed","4PL"), pch=c(19,-1), lty=c(-1,1), lwd=3,
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