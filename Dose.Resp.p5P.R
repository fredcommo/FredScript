Dose.Resp.p5P <- function(x, y, sdev, w.coef=1){

# x		: doses en log10, except 0
# y		: moyenne des réponses, en prop du Ctrl
# sdev	: sd des réponses.
# w.coef	: coefficient de poids

# Avec nlm
	# Fonction logistique 5P
		nlm.fit <- function(top, xmid, scal, s,  X){
			Y.fit <-(top)/(1+10^((xmid-X)*scal))^s
			return(Y.fit)
			}

	# Fonction sce (somme carré résidus) avec pondérations
		sce <- function(param, X, yobs, Weights) {
			top <- param[1]
			xmid <- param[2]
			scal <- param[3]
			s <- param[4]
			ytheo <- nlm.fit(top, xmid, scal, s, X)
			# residus <- yobs - ytheo
			# Weights <- sqrt(1/(residus^2))
			return(sum(Weights*(yobs - ytheo)^2))
			}


	# initialisation des valeurs
		z <- y/(max(y)*1.05)
		e = 1e-5
		weights <- ((z+e)*(1-z+e))^(1/w.coef)

		top.ini = max(y); max(y)
		xmid.ini = (max(x)+min(x))/2; xmid.ini
		scal.ini = 1/coef(lm(x~log(z/(1.005-z))))[2]
		scal.ini <- as.numeric(scal.ini); scal.ini
		s.ini = 1

	#rep(1, length(x))

		init <- c(top.ini, xmid.ini, scal.ini, s.ini)
		best<-nlm(f = sce, p = init, X = x, yobs = y, Weights=weights)						# calcul des paramètres

	# Récupération des paramètres
		best.top <- best$estimate[1]
		best.xmid<-best$estimate[2]
		best.scal <- best$estimate[3]
		best.s <- best$estimate[4]

	# Estimation des valeurs
		newX <- seq(min(x), max(x)*1.5, length=1000)						
		Y.best <- nlm.fit(best.top, best.xmid, best.scal, best.s, newX)

			# coordonnées du pt d'inflexion
				Xflex = best.xmid + (1/best.scal)*log10(best.s)
				Yflex = (best.top)*(best.s/(best.s+1))^best.s

			# coordonnées du pt rep = 0.5
				# Y50 = (max(y) + min(y))/2
				# Y50 = (best.top)/2
				Y50 = 0.5
				X50 = best.xmid - 1/best.scal*log10(((best.top)/(Y50))^(1/best.s)-1)
				# X50 = best.xmid - 1/best.scal*log10((1/(Y50))^(1/best.s)-1)								# top = 1, bottom = 0

			# pente au pt d'inflexion
				# B = (best.top)*log(10)*(best.scal)*(best.s/(best.s+1))^(best.s+1); B	# + best.bottom 							# 5P	 	pour d=1, (d/(d+1))^(d+1) = 1/4
				B = log(10)*(best.scal)*(best.s/(best.s+1))^(best.s+1); B								# top = 1, bottom = 0
				A = Yflex  - B*(Xflex); A													# 5P		pour d=1, (d/(d+1))^d = 1/2
				print(c(top = best.top,
						xmid = best.xmid,
						scal = best.scal,
						s = best.s,
						Xflex = Xflex, Yflex = Yflex,
						X50 = X50, D50 = 10^X50))

	# Représentations graphiques
		plot(y~x, pch = 20, cex = 1.5, col = "blue", xlim=range(min(newX), max(newX)),ylim = range(0,1), xlab="Log10(dose)", ylab="Response vs. Ctrl")
		lines(Y.best~newX, col = "orangered", lwd = 2)
		# lines(I(A+B*newX)~newX, col="lightblue", lty = 3, lwd = 1)
		abline(h = 0, lwd = 2)
		# points(Xflex, Yflex, pch = 19, col = "lightblue2", cex = 1.25)
		points(X50, Y50, pch = 8, col = "black", cex = 1.5)
		segments(x0 = X50, y0 = 0, x1 = X50, y1 = Y50, lty = 3, col="orangered")
		segments(x0 = min(newX)*2, y0 = Y50, x1 = X50, y1 = Y50, lty = 3, col="orangered")
		points(-A/B, 0, pch = 19, col = "red")

		# A compléter: sdev <- sd(as.data.frame(t(Resp[-1,])), na.rm = T)		

		for (i in 1:length(y)){
			pas <- (max(x)+min(x))/2*0.5
			segments(x0=x[i], y0=y[i]-sdev[i], x1=x[i], y1=y[i]+sdev[i], lty=3)
			segments(x0=x[i]-pas, y0=y[i]-sdev[i], x1=x[i]+pas, y1=y[i]-sdev[i], lty=3)
			segments(x0=x[i]-pas, y0=y[i]+sdev[i], x1=x[i]+pas, y1=y[i]+sdev[i], lty=3)
			}

		list(top = best.top, xmid = best.xmid, scal = best.scal, s = best.s, Xflex = Xflex, Yflex = Yflex, X50 = X50, D50 = 10^X50)

}