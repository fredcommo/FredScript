IC50.4P <- function(dose, rep, sd= NULL, W.coef=0.5, unit="nM", Title="", Plot.sd = TRUE, AddIC = TRUE, Add = FALSE, pCol = "grey50", lCol = "grey25",...){



	# Fonction logistique 4P
		nlm.fit.4P <- function(bottom, top, xmid, scal, X){
			Y.fit <-bottom+(top-bottom)/(1+10^((xmid-X)*scal))
			return(Y.fit)
			}

	# Fonction sce (somme carré résidus) avec pondérations
		sce.4P <- function(param, X, yobs, Weights, W.coef) {
			bottom <- param[1]
			top <- param[2]
			xmid <- param[3]
			scal <- param[4]
			ytheo <- nlm.fit.4P(bottom, top, xmid, scal, X)
			residus <- yobs - ytheo
			Weights <- (1/(residus^2))^(W.coef)
			return(sum(Weights*(yobs - ytheo)^2))
			}

	# initialisation des valeurs
		# add.first = FALSE			# Ajoute un point suppl. ?
		# if(y[1]<0.9){
			x.marge <- (max(dose)-min(dose))/(length(dose)-1)
			x <- c(min(dose)-x.marge, dose)
			y <- c(1, rep)
			sd <- c(0, sd)
			add.first = TRUE
		#	}

		w.coef = W.coef
		bottom.ini = miny = min(y, na.rm = T)	# ; min(y)
		top.ini = maxy = max(y, na.rm = T)		# ; max(y)
		maxx = max(x, na.rm = T)
		minx = min(x, na.rm = T)
		xmid.ini = (maxx + minx)/2; xmid.ini
		z <- (y - bottom.ini)/(top.ini - bottom.ini)	# z <- y/max(y, na.rm = T)
		z[z==0] <- 0.01
		z[z==1] <- 0.99
		scal.ini = coef(lm(x~log(z/(1-z))))[2]
		scal.ini <- as.numeric(scal.ini); scal.ini
		e = 1e-5
		# weights <- (z*(1-z))^(w.coef)
		# weights <- weights/sum(weights)
		weights <- rep(1, length(y))		

		init <- c(bottom=bottom.ini, top=top.ini, xmid=xmid.ini, scal=scal.ini)
		best<-nlm(f = sce.4P, p = init, X = x, yobs = y, Weights = weights, W.coef = w.coef)						# calcul des paramètres

	# Récupération des paramètres
		best.bottom <- best$estimate[1]
		best.top <- best$estimate[2]
		best.xmid<- best$estimate[3]
		best.scal <- best$estimate[4]

	# Estimation des valeurs
		newX <- seq(minx - x.marge, maxx + x.marge, length = 100)						
		Y.best <- nlm.fit.4P(best.bottom, best.top, best.xmid, best.scal, newX)

			# coordonnées du pt d'inflexion
				Xflex = best.xmid + (1/best.scal)*log10(1)
				Yflex = best.bottom + (best.top - best.bottom)*(1/(1+1))^1

			# coordonnées du pt rep = 0.5
				Y50 = 0.5
				X50 = best.xmid - 1/best.scal*log10(((best.top - best.bottom)/(Y50 - best.bottom))^(1/1)-1)

			# pente au pt d'inflexion
				B = best.bottom + (best.top - best.bottom)*log(10)*(best.scal)*(1/(1+1))^(1+1); B	# + best.bottom 							# 5P	 	pour d=1, (d/(d+1))^(d+1) = 1/4
				A = Y50  - B*(X50); A													# 5P		pour d=1, (d/(d+1))^d = 1/2

			# Intervalle IC X50
				ytheo <- nlm.fit.4P(best.bottom, best.top, best.xmid, best.scal, x)
				Q <- sum((y - ytheo)^2)
				S2 <- Q/(length(x)-2)
				s <- sqrt(S2)

				MC.X50 <- rep(0, 10000)

				for (i in 1:10000)
					MC.X50[i] <- best.xmid - 1/best.scal*log10(((best.top - best.bottom)/(Y50+rnorm(1,0,s) - best.bottom))^(1/1)-1)

				min.X50 <- quantile(MC.X50, probs=0.025, na.rm=T)
				min.D50 <- signif(10^min.X50, 2)
				max.X50 <- quantile(MC.X50, probs=0.975, na.rm=T)
				max.D50 <- signif(10^max.X50, 2)



	# Représentations graphiques
		if(!Add) plot(y[-1] ~ x[-1], pch = 20, cex = 2.5, col = pCol, xlim = range(minx, maxx),				# , ylim = range(0,max(maxy, Y.best)) , ylab = "Survival (% of control)"
					xlab = paste("Log10[", unit, "]", sep = ""),...)
		else points(y[-1] ~ x[-1], pch = 20, cex = 2.5, col = pCol, xlim = range(minx, maxx), ylim = range(0,max(maxy, Y.best)))

		# if (add.first) points(y[1] ~ x[1], pch = 20, cex = 2.5, col = "grey")

		if(AddIC){
			legend1 <- paste("IC50 :", round(10^X50, 2),unit)
			legend2 <- paste("[", min.D50, ",", max.D50, "]")
			# legend("bottomleft", legend = c(legend1, legend2), bty="n")
			text(c(x[1], x[1]), c(0.3, 0.2), pos=4, labels = c(legend1, legend2), font=4, cex=1.4, col="royalblue4")
			}

		if(Plot.sd){
			pas = x.marge/4
			for(s in 2:length(sd)){
				segments(x0 = x[s], x1 = x[s], y0 = y[s]-sd[s], y1 = y[s]+sd[s], lty = 3)
				segments(x0 = x[s]-pas, x1 = x[s]+pas, y0 = y[s]-sd[s], y1 = y[s]-sd[s], lty = 3)
				segments(x0 = x[s]-pas, x1 = x[s]+pas, y0 = y[s]+sd[s], y1 = y[s]+sd[s], lty = 3)
				}
			}

		lines(Y.best ~ newX, col = lCol, lwd = 3)
		Sub = "Weighted 4P logistic regr. (F. Commo)"
		if(W.coef==0) Sub = "Non weighted 4P logistic regr."
		title (main = Title, sub = Sub)

	# Output
		return(list(bottom = best.bottom, top = best.top,
						xmid = best.xmid, scal = best.scal,
						Xflex = Xflex, Yflex = Yflex, X50 = X50,
						D50 = round(10^X50, 2), Dmin = min.D50, Dmax = max.D50,
						xFit = newX, yFit = Y.best))

}


