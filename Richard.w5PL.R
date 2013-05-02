Richard.w5PL <- function(x, y, w = 0.25, Target = 0.5, Plot = T, add.line = T, add.points = T, add.intercept = T, pcol = "royalblue1", lcol = "navy", tan.col = "purple", Xlim = range(x), Ylim = range(y),...){

# x				: x-axis values
# y				: y-axis values
# w = 0.25			: weights coefficient
# Plot = F			: if results have to be visualized
# add.points = F		: add points on the plot
# add.intercept = T	: add the intercept point on plot
# pcol = "royalblue1"	: define points color	
# add.line = F		: add the tangente line
# lcol = "navy"		: define the regression color
# tan.col = "purple"	: define the tangente line color
# Xlim = range(x)		: the x-axis range
# Ylim = range(y)		: the y-axis range
# Title = ""		: to add a plot title


x <- as.numeric(x)
y <- as.numeric(y)

if(any(is.na(x), is.na(y))){
	na.index <- which(is.na(x) | is.na(y))
	y <- y[-na.index]
	x <- x[-na.index]
	}

# Fonction logistique 5PL
	Richard <- function(x, Fb, Fmax, b, c, d){
		y <- Fb + (Fmax-Fb)/(1 + exp(-(x-c)/b))^d
		return(y)
		}

# Fonction sce (somme carré résidus) avec pondérations
	sce.5P <- function(param, xobs, yobs, Weights, w) {
		Fb <- param[1]
		Fmax <- param[2]
		b <- param[3]
		c <- param[4]
		d <- param[5]
		ytheo <- Richard(xobs, Fb, Fmax, b, c, d)
		sq.res <- (yobs - ytheo)^2
		if(any(sq.res == 0)) sq.res[sq.res == 0] <- 1e-3
		Weights <- (1/sq.res)^w
		return(sum(Weights*sq.res))
		}

# Fonction sce (somme carré résidus) avec pondérations
	sce.5P.diag <- function(yobs, ytheo, w) {
		sq.res <- (yobs - ytheo)^2
		weights <- (1/sq.res)^w
		return(weights)
		}

# initialisation des parametres
Fb.ini = min(y, na.rm = T)
Fmax.ini = max(y, na.rm = T)	#*1.05
c.ini = (max(x, na.rm = T) + min(x, na.rm = T))/2
z <- (y - Fb.ini)/(Fmax.ini - Fb.ini)
z[z>=1] <- 0.95
z[z<=0] <- 0.05
b.ini = 1/coef(lm(log(z/(1-z))~x))[2]	# remplacer par b.ini = 1
# b.ini = 1
d.ini = sqrt(Fmax.ini)
init <- c(Fb.ini, Fmax.ini, b.ini, c.ini, d.ini)
weights <- rep(1, length(y))

# Estimation du modele
	best <- nlm(f = sce.5P, p = init, xobs = x, yobs = y, Weights = weights, w = w)

# Récupération des paramètres
	best.Fb <- best$estimate[1]
	best.Fmax <- best$estimate[2]
	best.b <- best$estimate[3]
	best.c <- best$estimate[4]
	best.d <- best$estimate[5]

# Diagnostic de régression
	yfit <- Richard(x, best.Fb, best.Fmax, best.b, best.c, best.d)
	weights <- sce.5P.diag(y, yfit, w)
	lm.test <- lm(yfit~y)	#, weights = weights)
	r.sq <- summary(lm.test)$adj.r.squared
	lm.slope <-coef(lm.test)[2]
	p.slope <- summary(lm.test)$coefficients[2,4]

# Estimation des valeurs pour graphique
	newx <- seq(min(x), max(x), length=100)						
	newy <- Richard(newx, best.Fb, best.Fmax, best.b, best.c, best.d)

# coordonnées du pt d'inflexion
	Xflex = best.c + best.b*log(best.d)
	Yflex = best.Fb + (best.Fmax-best.Fb)*(best.d/(1 + best.d))^(best.d)

# coordonnées du pt Target
	Y50 = min(newy) + (max(newy) - min(newy))*Target
	X50 = best.c - best.b*log(((best.Fmax-best.Fb)/(Y50 - best.Fb))^(1/best.d) - 1)

# pente au pt d'inflexion
	B = (best.Fmax-best.Fb)/best.b*(best.d/(1 + best.d))^(best.d + 1); B
	A = Yflex  - B*(Xflex); A

# pente finale
	x.ini <- x[1]
	x.end <- x[length(x)]
	y.ini <- Richard(x.ini, best.Fb, best.Fmax, best.b, best.c, best.d)
	y.end <- Richard(x.end, best.Fb, best.Fmax, best.b, best.c, best.d)
	Bf = (best.d*(best.Fmax-best.Fb)/best.b)*exp(-1/best.b*(x.end-best.c))*(1+exp(-1/best.b*(x.end-best.c)))^(-best.d-1)
	Af = y.end - Bf*x.end
 
# x.intercept
Cy0 = -A/B

# Représentations graphiques
	if(Plot){
	plot(y~x, col = pcol,...)
	lines(newy~newx, col = "navy", lwd = 4)
	if(add.points) points(y~x, col = pco,...l)
	if(add.line){
		lines(I(A+B*x)~x, col = tan.col, lwd = 1)
		abline(h = 0, lwd = 2)
		}
	if(add.intercept) points(-A/B, 0, pch = 19, col = "red")
	}

return(list(bottom = best.Fb, top = best.Fmax, xmid = best.c, scal = best.b, d = best.d,
			Xflex = Xflex, Yflex = Yflex, Xtarget = X50, Ytarget = Y50,
			slope = B, x.intercept = Cy0, Yini = y.ini, Yend = y.end, end.slope = Bf,
			lm.rsq = r.sq, lm.slope = lm.slope, p.value = p.slope, xfit = newx, yfit = newy))
}

