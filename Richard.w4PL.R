Richard.w4PL <- function(x, y, w = 0.25, Plot = F, add.points = F, add.intercept = T, pcol = "royalblue1", add.line = F, lcol = "navy", tan.col = "purple", Xlim = range(x), Ylim = range(y), Title = ""){

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

if(any(is.na(y))){
	na.index <- which(is.na(y))
	y <- y[-na.index]
	x <- x[-na.index]
	}

# Fonction logistique 5PL
	Richard <- function(x, Fb, Fmax, b, c){
		y <- Fb + Fmax/(1 + exp(-(x-c)/b))
		return(y)
		}

# Fonction sce (somme carré résidus) avec pondérations
	sce.4P <- function(param, xobs, yobs, w) {
		Fb <- param[1]
		Fmax <- param[2]
		b <- param[3]
		c <- param[4]
		ytheo <- Richard(xobs, Fb, Fmax, b, c)
		sq.res <- (yobs - ytheo)^2
		weights <- 1/sq.res^w
		return(sum(weights*sq.res))
		}

# Fonction sce (somme carré résidus) avec pondérations
	sce.5P.diag <- function(yobs, ytheo, w) {
		sq.res <- (yobs - ytheo)^2
		weights <- 1/sq.res^w
		return(weights)
		}

# initialisation des parametres
Fb.ini = min(y)
Fmax.ini = max(y)	#*1.05
c.ini = (max(x) + min(x))/2
z <- (y)/(Fmax.ini - y)
z <- abs(z)
if (any(abs(z)==Inf)) z[abs(z)==Inf] <- NA
b.ini = coef(lm(x~log(z)))[2]
# b.ini = 1					
init <- c(Fb.ini, Fmax.ini, b.ini, c.ini)

# Estimation du modele
	best<-nlm(f = sce.4P, p = init, xobs = x, yobs = y, w = w)

# Récupération des paramètres
	best.Fb <- best$estimate[1]
	best.Fmax <- best$estimate[2]
	best.b<-best$estimate[3]
	best.c <- best$estimate[4]

# Diagnostic de régression
	yfit <- Richard(x, best.Fb, best.Fmax, best.b, best.c)
	weights <- sce.5P.diag(y, yfit, w)
	lm.test <- lm(yfit~y)	#, weights = weights)
	r.sq <- summary(lm.test)$adj.r.squared
	lm.slope <-coef(lm.test)[2]
	p.slope <- summary(lm.test)$coefficients[2,4]

# Estimation des valeurs pour graphique
	newx <- seq(min(x), max(x), length=100)						
	newy <- Richard(newx, best.Fb, best.Fmax, best.b, best.c)

# coordonnées du pt d'inflexion
	best.d = 1
	Xflex = best.c + best.b*log(best.d)
	Yflex = best.Fb + best.Fmax*(best.d/(1 + best.d))^(best.d)

# coordonnées du pt rep = 0.5
	Y50 = (max(newy) + min(newy))/2
	X50 = best.c - best.b*log((best.Fmax/(Y50 - best.Fb))^(1/best.d) - 1)

# pente au pt d'inflexion
	B = best.Fmax/best.b*(best.d/(1 + best.d))^(best.d + 1); B
	A = Yflex  - B*(Xflex); A

# pente finale
	x.ini <- x[1]
	x.end <- x[length(x)]
	y.ini <- Richard(x.ini, best.Fb, best.Fmax, best.b, best.c)
	y.end <- Richard(x.end, best.Fb, best.Fmax, best.b, best.c)
	Bf = (best.d*best.Fmax/best.b)*exp(-1/best.b*(x.end-best.c))*(1+exp(-1/best.b*(x.end-best.c)))^(-best.d-1)
	Af = y.end - Bf*x.end
 
# x.intercept
Cy0 = -A/B

# Représentations graphiques
	if(Plot){
	plot(y~x, pch = 19, cex = 1.25, col = pcol, xlim = Xlim, ylim = Ylim, main = Title)
	lines(newy~newx, col = "navy", lwd = 2)
	lines(I(A+B*x)~x, col = tan.col, lwd = 1)
	abline(h = 0, lwd = 1, lty = 3)
	if(add.intercept) points(-A/B, 0, pch = 19, col = "red")
	}

	if(add.points) points(y~x, pch = 19, cex = 1.25, col = pcol)

	if(add.line){
	lines(newy~newx, col = lcol, lwd = 2)
	lines(I(A+B*x)~x, col = tan.col, lwd = 1)
	abline(h = 0, lwd = 2)
	if(add.intercept) points(-A/B, 0, pch = 19, col = "red")
	}

return(list(bottom = best.Fb, top = best.Fmax, xmid = best.c, scal = best.b,
			Xflex = Xflex, Yflex = Yflex, slope = B, x.intercept = Cy0, Yini = y.ini, Yend = y.end, end.slope = Bf,
			lm.rsq = r.sq, lm.slope = lm.slope, p.value = p.slope, xfit = newx, yfit = newy))
}

