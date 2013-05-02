Richard.5PL <- function(x, y, Plot = F){

	if(any(is.na(y))){
		na.index <- which(is.na(y))
		y <- y[-na.index]
		x <- x[-na.index]
		}

	# Fonction logistique 5PL
		Richard <- function(x, Fb, Fmax, b, c, d){
			y <- Fb + (Fmax - Fb)/(1 + exp(-(x-c)*b))^d
			return(y)
			}

	# Fonction sce (somme carré résidus) avec pondérations
		sce.5P <- function(param, xobs, yobs) {
			Fb <- param[1]
			Fmax <- param[2]
			b <- param[3]
			c <- param[4]
			d <- param[5]
			ytheo <- Richard(xobs, Fb, Fmax, b, c, d)
			return(sum((yobs - ytheo)^2))
			}

	Fb.ini = min(y); Fb.ini
	Fmax.ini = max(y); Fmax.ini
	c.ini = (max(x) + min(x))/2; c.ini
	z <- y/(Fmax.ini*1.05)
	b.ini = 1/coef(lm(x~log(z/(1-z))))[2]
	b.ini <- as.numeric(b.ini); b.ini
	d.ini = 1

	init <- c(Fb.ini, Fmax.ini, b.ini, c.ini, d.ini)

	# Estimation du modele
		best<-nlm(f = sce.5P, p = init, xobs = x, yobs = y)

	# Récupération des paramètres
		best.Fb <- best$estimate[1]
		best.Fmax <- best$estimate[2]
		best.b<-best$estimate[3]
		best.c <- best$estimate[4]
		best.d <- best$estimate[5]

	# Estimation des valeurs
		newx <- seq(min(x), max(x), length=100)						
		yfit <- Richard(newx, best.Fb, best.Fmax, best.b, best.c, best.d)

	# coordonnées du pt d'inflexion
		Xflex = best.c + 1/best.b*log(best.d)
		Yflex = best.Fb + (best.Fmax - best.Fb)*(best.d/(1 + best.d))^(best.d)

	# coordonnées du pt rep = 0.5
		Y50 = (max(yfit) + min(yfit))/2
		X50 = best.c - 1/best.b*log(((best.Fmax - best.Fb)/(Y50 - best.Fb))^(1/best.d) - 1)

	# pente au pt d'inflexion
	B = (best.Fmax - best.Fb)*best.b*(best.d/(1 + best.d))^(best.d + 1); B
	A = Yflex  - B*(Xflex); A

	# x.intercept
	Cy0 = -A/B

	# Représentations graphiques
		if(Plot){
			plot(y~x, pch = 1, cex = 1.25, col = "royalblue1", xlim = range(newx), ylim = range(yfit))
			lines(yfit~newx, col = "navy", lwd = 2)
			lines(I(A+B*x)~x, col="purple", lwd = 1)
			abline(h = 0, lwd = 2)
			points(-A/B, 0, pch = 19, col = "red")
			}

	return(list(bottom = best.Fb, top = best.Fmax, xmid = best.c, scal = best.b, d = best.d,
						Xflex = Xflex, Yflex = Yflex, slope = B, x.intercept = Cy0))

}
