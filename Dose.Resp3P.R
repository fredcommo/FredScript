Dose.Resp3P <- function(Dose, Resp, Prop=FALSE, Cor=FALSE, Blank=0, T0=0, wi=T, w.coef=2, add=F, print.legend=TRUE, ExcludeMinMax=FALSE, unit="µM"){

# Dose : vecteur doses
# Resp : data frame des valeurs
# Prop : si TRUE, calcule en prop max/min et renvoie l'EC50
# Cor : si TRUE, les valeurs sont les valeurs corrigées et en prop du ctrl.

if (ExcludeMinMax){
	i.min <- apply(Resp, 1, which.min)
	i.max <- apply(Resp, 1, which.max)
		for (i in 1:nrow(Resp)) Resp[i,c(i.min[i], i.max[i])] <- NA
	}

	ctrl.index <- which(Dose==0)
	x <- Dose[-ctrl.index]
	x <- log10(x)

	y <- mean(as.data.frame(t(Resp[-ctrl.index,])), na.rm = T)
	sdev <- sd(as.data.frame(t(Resp[-ctrl.index,])), na.rm = T)

	if(!Cor){
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
		}

	# Avec nlm
		# Fonction logistique 3P
			nlm.fit <- function(top, xmid, scal, X){
				Y.fit <- top/(1+10^((xmid-X)*scal))
				return(Y.fit)
				}

		# Fonction sce (somme carré résidus) avec pondérations
			sce <- function(param, X, yobs, Weights) {
				top <- param[1]
				xmid <- param[2]
				scal <- param[3]
				ytheo <- nlm.fit(top, xmid, scal, X)
				# residus <- (yobs - ytheo)^2
				# Weights <- sqrt(1/(residus))
				return(sum(Weights*(yobs - ytheo)^2))
				}
	
		# initialisation des valeurs
			top.ini = max(y, na.rm=T)
			xmid.ini = (max(x)+min(x))/2
			z <- y	#/(top.ini*1.05)						# y sont les prop vs. Ctrl => z sont les prop vs. prop max
	
			if (Prop){
				min.y <- min(y, na.rm=T)*0.995
				max.y <- max(y, na.rm=T)*1.005
				z <- (y-min.y)/(max.y-min.y)			# si y = prop
				}
	
			scal.ini = 1/coef(lm(x~log(z/(1.005-z))))[2]
			scal.ini <- as.numeric(scal.ini)
			weights <- rep(1, length(z))
		
			if (wi) weights <- (z*(1-z))^(1/w.coef)
	
			init <- c(top=top.ini, xmid=xmid.ini, scal=scal.ini)
			best<-nlm(f = sce, p = init, X = x, yobs = z, Weights=weights)						# calcul des paramètres
	
		# Récupération des paramètres
			best.top <- best$estimate[1]
			best.xmid<- best$estimate[2]
			best.scal <- best$estimate[3]	

		# Estimation des valeurs
			newX <- seq(min(x), max(x)*1.5, length=100)						
			Y.best <- nlm.fit(best.top, best.xmid, best.scal, newX)
	
				# coordonnées du pt d'inflexion
					Xflex = best.xmid + (1/best.scal)*log10(1)
					Yflex = best.top*(1/(1+1))^1

				# coordonnées du pt rep = 0.5
					# Y50 = (max(y) + min(y))/2
					# Y50 = (best.bottom + best.top)/2
					Y50 = 0.5
					X50 = best.xmid - 1/best.scal*log10((best.top/Y50)^(1/1)-1)
					# X50 = best.xmid - 1/best.scal*log10((1/(Y50))^(1/1)-1)								# top = 1, bottom = 0

				# pente au pt d'inflexion
					B = best.top*log(10)*(best.scal)*(1/(1+1))^(1+1); B	# + best.bottom 							# 5P	 	pour d=1, (d/(d+1))^(d+1) = 1/4
					# B = log(10)*(best.scal)*(1/(1+1))^(1+1); B								# top = 1, bottom = 0
					A = Y50  - B*(X50); A													# 5P		pour d=1, (d/(d+1))^d = 1/2
					print(c(top = best.top,
							xmid = best.xmid,
							scal = best.scal,
							Xflex = Xflex, Yflex = Yflex,
							X50 = X50, Y50 = Y50))

				# Intervalle IC X50
					ytheo <- nlm.fit(best.top, best.xmid, best.scal, x)
					Q <- sum((y - ytheo)^2)
					S2 <- Q/(length(x)-2)
					sd <- sqrt(S2)
	
					MC.X50 <- rep(0, 10000)
	
					for (i in 1:10000)
						MC.X50[i] <- best.xmid - 1/best.scal*log10((best.top/(Y50+rnorm(1,0,2*sd)))^(1/1)-1)

					q.X50 <- quantile(MC.X50, probs=c(0.025, 0.975), na.rm=T)



		# Représentations graphiques
			pt.col= lin.col = "orangered"
			if (add){ 
				pt.col = lin.col = "blue"
				points(z~x, pch = 20, cex = 1.25, col = pt.col)
				}

			if (!add)
			plot(z~x, pch = 20, cex = 1.25, col = pt.col, 
				xlim=range(min(newX),max(newX)), ylim = range(0, 1),
				xlab="Log10(drug)", ylab="Survival vs. control")


			for (i in 1:length(y)){
				pas <- (max(x)+min(x))/2*0.5
				segments(x0=x[i], y0=y[i]-sdev[i], x1=x[i], y1=y[i]+sdev[i], lty=3)
				segments(x0=x[i]-pas, y0=y[i]-sdev[i], x1=x[i]+pas, y1=y[i]-sdev[i], lty=3)
				segments(x0=x[i]-pas, y0=y[i]+sdev[i], x1=x[i]+pas, y1=y[i]+sdev[i], lty=3)
				}

			lines(Y.best~newX, col = lin.col, lwd = 2)
			#segments(x0 = X50, y0 = 0, x1 = X50, y1 = Y50, lty = 3)
			#segments(x0 = 0, y0 = Y50, x1 = X50, y1 = Y50, lty = 3)

			res <- list(D50=10^X50, Dmin=10^q.X50[1], Dmax=10^q.X50[2])
			results <- paste("IC50 =", signif(res$D50,3), unit, "[", signif(res$Dmin,3), ",", signif(res$Dmax,3), "]")

			if (Prop)
				results <- paste("EC50 =", signif(res$D50,3), unit, "[", signif(res$Dmin,3), ",", signif(res$Dmax,3), "]")

			if(print.legend) legend("topright", legend=results, bty="n")
	
			print(res)
}		
		# Ajouter l'intervalle de confiance sur la courbe
