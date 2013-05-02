Quadratic <- function(x, y, Plot=TRUE, xlab="", ylab="", title=""){


# Calcul de l'alignement
	nlm.fit<-function(bottom, top, s, C, X){
		Y.fit <- C + bottom*X + (top - bottom)*(X^s)
		return(Y.fit)
		}
	
# Fonction sce (somme carré résidus)
		sce <- function(param, X, yobs, weights) {
			bottom <- param[1]
			top <- param[2]
			s <- param[3]
			C <- param[4]
			ytheo <- nlm.fit(bottom, top, s, C, X)
			residus <- (yobs - ytheo)^2
			weights <- (1/residus)^(1/4)
			# weights <- weights/mean(weights)
			return(sum(weights*(yobs - ytheo)^2))
			}

	# parametres initiaux
	
		if(is.null(weights)) weights <- rep(1, length(Y))

		init.bottom = min(Y, na.rm=T)
		init.top = max(Y, na.rm=T)
		init.s = 1
		init.C = init.bottom
	
		# {
		# if (y[1]<y[length(y)])	{
		#	min.y <- min(Y)
		#	max.y <- max(Y)
		#	}
		# else	{
		#	min.y <- max(y)
		#	max.y <- min(y)
		#	}
		# }

		# p <- (Y - min.y)/(max.y - min.y); p[order(X)]
		# weights <- (1-p^(1-p))^(2)
		# weights <- rep(1, length(Y))

	
	# Y.init<-nlm.fit(init.bottom, init.top, init.k, init.X.offset, X)

		init=c(bottom=init.bottom, top=init.top, s=init.s, C=init.C)	# initialisation des paramètres
		best<-nlm(f = sce, p = init, X = X, yobs = Y, weights=weights)	# calcul des paramètres
	
		# Récupération des paramètres
			best.bottom <- best$estimate[1]
			best.top <- best$estimate[2]
			best.s <- best$estimate[3]
			best.C <- best$estimate[4]

		SStot <- length(X)*var(X)
		SSres <- sce(c(best.bottom, best.top, best.s, best.C), X, Y, weights=weights)
		r2 = (1 - SSres/SStot)

		newX <- seq(min(X), max(X), length=100)						
		Y.best <- nlm.fit(best.bottom, best.top, best.s, best.C, newX)
		
		slope = best.bottom + best.s*(best.top - best.bottom)*(newX[1])^(best.s-1); slope
		intercept = Y.best[1] - slope*newX[1]; intercept	

		if(Plot){
			y.range <- c(min(c(Y, Y.best)), max(c(Y, Y.best)))
			plot(Y ~ X, pch=19, col="blue", main="Quadratic Regression", ylim=y.range, sub=title, xlab=xlab, ylab=ylab)
			lines(I(intercept + slope*newX)~newX)
			lines(Y.best~newX, col="orangered", lwd=2)
			}

		test.lm <- lm.rob(nlm.fit(best.bottom, best.top, best.s, best.C, X), Y, Plot=F)
		
		list(Slope = slope, Intercept = intercept, Fit = test.lm, p.value=summary(test.lm)$coefficients[2, 4])
}