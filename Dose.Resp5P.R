Dose.Resp5P <- function(x, y, sdev, w.coef=1){

# Avec nlm
	# Fonction logistique 5P
		nlm.fit <- function(bottom, top, xmid, scal, s,  X){
			Y.fit <-bottom+(top-bottom)/(1+10^((xmid-X)*scal))^s
			return(Y.fit)
			}

	# Fonction sce (somme carré résidus) avec pondérations
		sce <- function(param, X, yobs, Weights) {
			bottom <- param[1]
			top <- param[2]
			xmid <- param[3]
			scal <- param[4]
			s <- param[5]
			ytheo <- nlm.fit(bottom, top, xmid, scal, s, X)
			# residus <- yobs - ytheo
			# Weights <- sqrt(1/(residus^2))
			return(sum(Weights*(yobs - ytheo)^2))
			}

	# initialisation des valeurs
		bottom.ini = min(y); min(y)
		top.ini = max(y); max(y)
		xmid.ini = (max(x)+min(x))/2; xmid.ini
		z <- y/(max(y)*1.05)
		# z <- y			# si y = prop
		scal.ini = 1/coef(lm(x~log(z/(1.005-z))))[2]
		scal.ini <- as.numeric(scal.ini); scal.ini
		s.ini = 1
		weights <- (z*(1-z))^(1/w.coef)

		init <- c(bottom.ini, top.ini, xmid.ini, scal.ini, s.ini)
		best<-nlm(f = sce, p = init, X = x, yobs = z, Weights=weights)						# calcul des paramètres

	# Récupération des paramètres
		best.bottom <- best$estimate[1]
		best.top <- best$estimate[2]
		best.xmid<-best$estimate[3]
		best.scal <- best$estimate[4]
		best.s <- best$estimate[5]

	# Estimation des valeurs
		newX <- seq(min(x), max(x)*1.5, length=1000)						
		Y.best <- nlm.fit(best.bottom, best.top, best.xmid, best.scal, best.s, newX)

			# coordonnées du pt d'inflexion
				Xflex = best.xmid + (1/best.scal)*log10(best.s)
				Yflex = best.bottom + (best.top - best.bottom)*(best.s/(best.s+1))^best.s

			# coordonnées du pt rep = 0.5
				# Y50 = (max(y) + min(y))/2
				# Y50 = (best.bottom + best.top)/2
				Y50 = 0.5
				X50 = best.xmid - 1/best.scal*log10(((best.top - best.bottom)/(Y50 - best.bottom))^(1/best.s)-1)
				# X50 = best.xmid - 1/best.scal*log10((1/(Y50))^(1/best.s)-1)								# top = 1, bottom = 0

			# pente au pt d'inflexion
				# B = best.bottom + (best.top - best.bottom)*log(10)*(best.scal)*(best.s/(best.s+1))^(best.s+1); B	# + best.bottom 							# 5P	 	pour d=1, (d/(d+1))^(d+1) = 1/4
				B = log(10)*(best.scal)*(best.s/(best.s+1))^(best.s+1); B								# top = 1, bottom = 0
				A = Yflex  - B*(Xflex); A													# 5P		pour d=1, (d/(d+1))^d = 1/2
				print(c(bottom = best.bottom, 
						top = best.top,
						xmid = best.xmid,
						scal = best.scal,
						s = best.s,
						Xflex = Xflex, Yflex = Yflex,
						X50 = X50, Y50 = Y50))

	# Représentations graphiques
		plot(z~x, pch = 20, cex = 0.75, col = "blue", ylim = range(0,1))
		lines(Y.best~newX, col = "orangered", lwd = 2)
		lines(I(A+B*x)~x, col="darkblue", lty = 3, lwd = 2)
		abline(h = 0, lwd = 2)
		points(Xflex, Yflex, pch = 19, col = "lightblue2", cex = 1.25)
		points(X50, Y50, pch = 8, col = "black", cex = 1.5)
		segments(x0 = X50, y0 = 0, x1 = X50, y1 = Y50, lty = 3)
		segments(x0 = 0, y0 = Y50, x1 = X50, y1 = Y50, lty = 3)
		points(-A/B, 0, pch = 19, col = "red")

}