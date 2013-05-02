IC50.5P.v2 <- function(x, y, sd = NULL, Weights = NULL, unit="nM", Title="", Plot.sd = TRUE, AddIC = TRUE, Add = FALSE, ptCol="indianred3", lCol="royalblue4"){

# Weighted by residuals

# Avec nlm
	# Fonction logistique 5P
		nlm.fit.5P <- function(bottom, top, xmid, scal, s,  X){
			Y.fit <-bottom+(top-bottom)/(1+10^((xmid-X)*scal))^s
			return(Y.fit)
			}

	# Fonction sce (somme carré résidus) avec pondérations
		sce.5P <- function(param, X, yobs, Weights) {
			bottom <- param[1]
			top <- param[2]
			xmid <- param[3]
			scal <- param[4]
			s <- param[5]
			ytheo <- nlm.fit.5P(bottom, top, xmid, scal, s, X)
			residus <- yobs - ytheo
			return(sum((1/Weights^2)*(yobs - ytheo)^2))
			}

	# initialisation des valeurs
		# add.first = FALSE			# Ajoute un point suppl. ?
		# if(y[1]<0.9){
			x.marge <- (max(x)-min(x))/(length(x)-1)
			x <- c(min(x)-x.marge, x)
			y <- c(1, y)
			sd <- c(0, sd)
			add.first = TRUE
			{
			if (is.null(Weights)) weights <- rep(1, length(y)+1)
			else weights <- c(1, Weights)
			}
	#	}


		bottom.ini = min(y); min(y)
		top.ini = max(y); max(y)
		xmid.ini = (max(x)+min(x))/2; xmid.ini
		z <- y/max(y)
		z[z==0] <- 0.01
		z[z==1] <- 0.99

		scal.ini = coef(lm(x~log(z/(1-z))))[2]
		scal.ini <- as.numeric(scal.ini); scal.ini
		s.ini = 1
		# weights <- (z*(1-z))^(w.coef)
		# weights <- weights/sum(weights)
		

		init <- c(bottom.ini, top.ini, xmid.ini, scal.ini, s.ini)
		best<-nlm(f = sce.5P, p = init, X = x, yobs = y, Weights=weights)						# calcul des paramètres

	# Récupération des paramètres
		best.bottom <- best$estimate[1]
		best.top <- best$estimate[2]
		best.xmid<-best$estimate[3]
		best.scal <- best$estimate[4]
		best.s <- best$estimate[5]

	# Estimation des valeurs
		newX <- seq(min(x)*0.75, max(x)*1.25, length=100)						
		Y.best <- nlm.fit.5P(best.bottom, best.top, best.xmid, best.scal, best.s, newX)

			# coordonnées du pt d'inflexion
				Xflex = best.xmid + (1/best.scal)*log10(best.s)
				Yflex = best.bottom + (best.top - best.bottom)*(best.s/(best.s+1))^best.s

			# coordonnées du pt rep = 0.5
				Y50 = 0.5
				X50 = best.xmid - 1/best.scal*log10(((best.top - best.bottom)/(Y50 - best.bottom))^(1/best.s)-1)

			# pente au pt d'inflexion
				# B = best.bottom + (best.top - best.bottom)*log(10)*(best.scal)*(best.s/(best.s+1))^(best.s+1); B	# + best.bottom 							# 5P	 	pour d=1, (d/(d+1))^(d+1) = 1/4
				B = log(10)*(best.scal)*(best.s/(best.s+1))^(best.s+1); B								# top = 1, bottom = 0
				A = Yflex  - B*(Xflex); A													# 5P		pour d=1, (d/(d+1))^d = 1/2

			# Intervalle IC X50
				ytheo <- nlm.fit.5P(best.bottom, best.top, best.xmid, best.scal, best.s,  x)
				Q <- sum((y - ytheo)^2)
				S2 <- Q/(length(x)-2)
				s <- sqrt(S2)

				MC.X50 <- rep(0, 10000)

				for (i in 1:10000){
					yi <- Y50 + rnorm(1, 0, s)
					MC.X50[i] <- best.xmid - 1/best.scal*log10(((yi-best.bottom)/(best.top-best.bottom))^(-1/best.s) - 1)
					}
				min.X50 <- quantile(MC.X50, probs=0.025, na.rm=T)
				min.D50 <- signif(10^min.X50, 2)
				max.X50 <- quantile(MC.X50, probs=0.975, na.rm=T)
				max.D50 <- signif(10^max.X50, 2)

	# Représentations graphiques
		if(!Add) plot(y[-1] ~ x[-1], pch = 20, cex = 2.5, col = ptCol, xlim = range(min(x), max(x)), ylim = range(0,max(y, Y.best)), xlab="Log10(Drug)", ylab="Survival (% of control)")
		else points(y[-1] ~ x[-1], pch = 20, cex = 2.5, col = ptCol, xlim = range(min(x), max(x)), ylim = range(0,max(y, Y.best)))

		# if (add.first) points(y[1] ~ x[1], pch = 20, cex = 2.5, col = "grey")

		if(AddIC){
			legend1 <- paste("IC50 :", round(10^X50, 2),unit)
			legend2 <- paste("[", min.D50, ",", max.D50, "]")
			# legend("bottomleft", legend = c(legend1, legend2), bty="n")
			text(c(x[1], x[1]), c(0.3, 0.2), pos=4, labels = c(legend1, legend2), font=4, cex=1.4, col="royalblue4")
			}

		if(Plot.sd){
			pas = x.marge/3
			for(s in 2:length(sd)){
				segments(x0 = x[s], x1 = x[s], y0 = y[s]-sd[s], y1 = y[s]+sd[s], lty = 3)
				segments(x0 = x[s]-pas, x1 = x[s]+pas, y0 = y[s]-sd[s], y1 = y[s]-sd[s], lty = 3)
				segments(x0 = x[s]-pas, x1 = x[s]+pas, y0 = y[s]+sd[s], y1 = y[s]+sd[s], lty = 3)
				}
			}

		lines(Y.best ~ newX, col = lCol, lwd = 3)
		Sub = "Weighted 5P logistic regr."
		if(is.null(Weights)) Sub = "Non weighted 5P logistic regr."
		title (main = Title, sub = Sub)

		# Output
		return(list(top = best.top, xmid = best.xmid, scal = best.scal, s = best.s,
				Xflex = Xflex, Yflex = Yflex, X50 = X50, 
				D50 = 10^X50, Dmin = min.D50, Dmax = max.D50,
				xFit = newX, yFit = Y.best))

}
