HillFit4 <- function (X, Y, test.rlm=TRUE, Slope=TRUE, Adjust=c("none","Scale","Prop"), xlab="", ylab="", weights=NULL, title=NULL, print=TRUE){
	
	# Ajustement des valeurs y
	ajust <- function(y, type=c("none","Scale","Prop")){
		type <- match.arg(type)
		switch(type, none=(y=y), Scale=scale(y), Prop=(y=y/(max(y, na.rm=T)*1.01)))
		}

	if (length(Adjust)==3) Adjust="none"
	Y <- ajust(Y, type=Adjust)

	# Calcul de l'alignement
		nlm.fit<-function(bottom,top,slope,X50,X){
			Y.fit<-bottom+(top-bottom)/(1+10^((X50-X)*slope))
			return(Y.fit)
			}

	# Fonction sce (somme carré résidus)
		sce <- function(param, X, yobs, weights) {
			bottom<-param[1]
			top<-param[2]
			slope<-param[3]
			X50<-param[4]
			ytheo <- nlm.fit(bottom, top, slope, X50, X)
			residus <- yobs - ytheo
			weights <- sqrt(1/(residus^2))
			return(sum(weights*(yobs - ytheo)^2))
			}

	# parametres initiaux
	
		if(is.null(weights)) weights <- rep(1, length(Y))
		lin.mod <- rlm(Y ~ X)
		init.slope = coef(lin.mod)[2]
		init.bottom = min(Y, na.rm=T)
		init.top = max(Y, na.rm=T)
		#if (init.slope<0){
		#	init.bottom = max(Y, na.rm=T)
		#	init.top = min(Y, na.rm=T)
		#}
		init.X50 = (max(X, na.rm = T)+min(X, na.rm = T))/2

	plot(Y~X, pch=19, col="blue", main="4P Logistic Regression", sub=title, xlab=xlab, ylab=ylab)
	
	# Y.init<-nlm.fit(init.bottom, init.top, init.slope, init.X50, X)

	cond = TRUE
	if (test.rlm) cond = abs(coef(lin.mod)[2])>0.005

	{
	if (cond) {

		init=c(bottom=init.bottom, top=init.top, slope=init.slope, X50=init.X50)	# initialisation des paramètres
		best<-nlm(f = sce, p = init, X = X, yobs = Y, weights=weights)						# calcul des paramètres
	
		# Récupération des paramètres
			best.bottom <- best$estimate[1]
			best.top <- best$estimate[2]
			Result <- best.slope<-best$estimate[3]
			best.X50 <- best$estimate[4]

		SStot <- length(X)*var(X)
		SSres <- sce(c(best.bottom, best.top, best.slope, best.X50), X, Y, weights=weights)
		r2 = (1-SSres/SStot)

		newX <- seq(min(X), max(X), length=100)						
		Y.best <- nlm.fit(best.bottom, best.top, best.slope, best.X50, newX)
		lines(Y.best~newX, col="orangered", lwd=2)
		
		if (Slope){
			Y50 <- (best.bottom + best.top)/2	#nlm.fit(best.bottom, best.top, best.slope, best.X50, best.X50)
			B = (best.top-best.bottom)*log(10)/(best.slope*4)
			A = (best.top + best.bottom)/2 - B*best.X50
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