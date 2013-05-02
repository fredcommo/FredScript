Dose.Resp4P.v2 <- function(Dose, Resp, Prop=FALSE, Blank=0, T0=0, wi=T, w.coef=1, add=F, print.legend=TRUE, ExcludeMinMax=FALSE, unit="µM"){

if (ExcludeMinMax){
	i.min <- apply(Resp, 1, which.min)
	i.max <- apply(Resp, 1, which.max)
		for (i in 1:nrow(Resp)) Resp[i,c(i.min[i], i.max[i])] <- NA
	}

ctrl.index <- which(Dose==0)
x <- Dose[-ctrl.index]
x <- log10(x)


	correct.resp <- Resp[-ctrl.index,] - T0
	ctrl <- mean(as.numeric(Resp[ctrl.index,]-T0), na.rm = T)

	y <- mean(as.data.frame(t(correct.resp/ctrl)), na.rm = T)
	sdev <- sd(as.data.frame(t(correct.resp/ctrl)), na.rm = T)

	whichNeg <- which(y<0)

	if(length(whichNeg)>0){	
		y <- y[-whichNeg]
		x <- x[-whichNeg]
		sdev <- sdev[-whichNeg]
		}

nrep <- ncol(Resp)

dose <- rep(x, nrep)
rep <- c()
for (i in 1:nrep) rep <- c(rep, Resp[,i])

# Avec nlm
	# Fonction logistique 4P
		nlm.fit <- function(bottom, top, xmid, scal, X){
			Y.fit <-bottom+(top-bottom)/(1+10^((xmid-X)*scal))
			return(Y.fit)
			}

	# Fonction sce (somme carré résidus) avec pondérations
		sce <- function(param, X, yobs, Weights) {
			bottom <- param[1]
			top <- param[2]
			xmid <- param[3]
			scal <- param[4]
			ytheo <- nlm.fit(bottom, top, xmid, scal, X)
			# residus <- (yobs - ytheo)^2
			# Weights <- sqrt(1/(residus))
			return(sum(Weights*(yobs - ytheo)^2))
			}

	# initialisation des valeurs
		bottom.ini = min(rep, na.rm=T); bottom.ini
		top.ini = max(rep, na.rm=T); top.ini
		xmid.ini = (max(dose)+min(dose))/2; xmid.ini
		z <- rep/(max(rep, na.rm=T)*1.05)

		if (Prop){
			min.rep <- min(rep, na.rm=T)*0.995
			max.rep <- max(rep, na.rm=T)*1.005
			z <- (rep-min.rep)/(max.rep - min.rep)			# si y = prop
			}

		scal.ini = 1/coef(lm(dose~log(z/(1.005-z))))[2]
		scal.ini <- as.numeric(scal.ini); scal.ini

		weights <- rep(1, length(z))
		if (wi) weights <- (z*(1-z))^(1/w.coef)		#rep(1, length(x))

		init <- c(bottom=bottom.ini, top=top.ini, xmid=xmid.ini, scal=scal.ini)
		best<-nlm(f = sce, p = init, X = dose, yobs = rep, Weights=weights)						# calcul des paramètres

	# Récupération des paramètres
		best.bottom <- best$estimate[1]
		best.top <- best$estimate[2]
		best.xmid<- best$estimate[3]
		best.scal <- best$estimate[4]

	# Estimation des valeurs
		newX <- seq(min(dose), max(dose)*1.5, length=1000)						
		Y.best <- nlm.fit(best.bottom, best.top, best.xmid, best.scal, newX)

			# coordonnées du pt d'inflexion
				Xflex = best.xmid + (1/best.scal)*log10(1)
				Yflex = best.bottom + (best.top - best.bottom)*(1/(1+1))^1

			# coordonnées du pt rep = 0.5
				# Y50 = (max(y) + min(y))/2
				# Y50 = (best.bottom + best.top)/2
				Y50 = 0.5
				X50 = best.xmid - 1/best.scal*log10(((best.top - best.bottom)/(Y50 - best.bottom))^(1/1)-1)
				# X50 = best.xmid - 1/best.scal*log10((1/(Y50))^(1/1)-1)								# top = 1, bottom = 0

			# pente au pt d'inflexion
				B = best.bottom + (best.top - best.bottom)*log(10)*(best.scal)*(1/(1+1))^(1+1); B	# + best.bottom 							# 5P	 	pour d=1, (d/(d+1))^(d+1) = 1/4
				# B = log(10)*(best.scal)*(1/(1+1))^(1+1); B								# top = 1, bottom = 0
				A = Y50  - B*(X50); A													# 5P		pour d=1, (d/(d+1))^d = 1/2
				print(c(bottom = best.bottom, 
						top = best.top,
						xmid = best.xmid,
						scal = best.scal,
						Xflex = Xflex, Yflex = Yflex,
						X50 = X50, Y50 = Y50))

			# Intervalle IC X50
				ytheo <- nlm.fit(best.bottom, best.top, best.xmid, best.scal, x)
				Q <- sum((y - ytheo)^2)
				S2 <- Q/(length(x)-2)
				sd <- sqrt(S2)

				MC.X50 <- rep(0, 10000)

				for (i in 1:10000)
					MC.X50[i] <- best.xmid - 1/best.scal*log10(((best.top - best.bottom)/(Y50+rnorm(1,0,sd) - best.bottom))^(1/1)-1)

				q.X50 <- quantile(MC.X50, probs=c(0.025, 0.975), na.rm=T)
				q.X50; 10^(q.X50)



	# Représentations graphiques
		plot(rep~dose, pch = 20, cex = 1.25, col = "blue", xlim=range(min(newX), max(newX)), ylim = range(0, 1))

		for (i in 1:length(y)){
			points(y ~ x, pch=19, cex=1.25, col="red")
			pas <- (max(x)+min(x))/2*0.05
			segments(x0=x[i], y0=z[i]-sdev[i], x1=x[i], y1=z[i]+sdev[i], lty=3)
			segments(x0=x[i]-pas, y0=z[i]-sdev[i], x1=x[i]+pas, y1=z[i]-sdev[i], lty=3)
			segments(x0=x[i]-pas, y0=z[i]+sdev[i], x1=x[i]+pas, y1=z[i]+sdev[i], lty=3)
			}

		lines(Y.best~newX, col = "orangered", lwd = 2)
		segments(x0 = X50, y0 = 0, x1 = X50, y1 = Y50, lty = 3)
		segments(x0 = 0, y0 = Y50, x1 = X50, y1 = Y50, lty = 3)

		print(data.frame(D50=10^X50, Dmin=10^q.X50[1], Dmax=10^q.X50[2]))

}


# Ajouter l'intervalle de confiance sur la courbe