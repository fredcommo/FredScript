myPairs <- function(data, smooth.fun = c("LM", "RLM", 'lmRob', "SMOOTH"), cor.method = c("spearman", "pearson"), Pcol = "royalblue3", Lcol = "red", CexCor = 3.5, subTitle = TRUE,...){

# data: a data set (or matrix) of observations by row, and variates by columns.
# smooth.fun: what smoothing function to use. "LM" performs and add a standard linear regression (default), "RLM" performs a robust linear regression, "SMOOTH" use a kernel smoothing.
# cor.method: what method to compute the correlation coefficients. Default use 'Spearman'
# Pcol: the points color.
# Lcol: the line color.
# CexCor: a parameter to adjust the size of the text in the upper panel.
# Title: if TRUE (default) a subtitle is added indicating the smoothing and correlation methods used.
#...: other optional graphic parameters.

require(MASS)

mylm <- function(x, y, Pcol, Lcol,...){
	lm.test <- lm(y~x)
	points(y~x, col = Pcol,...)
	abline(lm.test, col = Lcol,...)
	}

myrlm <- function(x, y, Pcol, Lcol,...){
	rlm.test <- rlm(y~x)
	points(y~x, col = Pcol,...)
	abline(rlm.test, col = Lcol,...)
	}

lmRob <- function (x, y, Pcol, Lcol, Eps = 1e-5, Coef = 5, ...) {
		na.index <- which(is.na(x) | x==Inf | x==(-Inf) | is.na(y) | y==Inf | y==(-Inf))
		if(length(na.index)>0) {x <- x[-na.index]; y <- y[-na.index]}
		epsilon = Eps
		lm.test <- lm(y ~ x)					#; abline(lm1)
		s <- summary(lm.test)$sigma
		s.init = Inf
		j = 0							# compteur iterations
		while(abs(s.init - s) > 1e-10 & j<1000){
			s.init <- s
			# Ponderations
				r <- resid(lm.test)
				m = median (r, na.rm = T)
				s = median (abs (r - m), na.rm = T)
				u = (r - m) / ((Coef * s) + epsilon)
				W = rep (0, length(x))
				i = abs(u) <= 1
				w <- rep(0, length(x))
				w[i] = ((1 - u^2)^2)[i]		
			# regression ponderee
			lm.test <- lm(y ~ x, weights = w)			#; abline(lm1)
			# nouveaux residus
			s <- summary(lm.test)$sigma
			j = j+1
			}
	points(y~x, col = Pcol,...)
	abline(lm.test, col = Lcol,...)
}		

panel.cor <- function(x, y, digits = 3, prefix="", mycor, cex.cor = CexCor){
		usr <- par("usr"); on.exit(par(usr))
		par(usr = c(0, 1, 0, 1))
		r <- cor(x, y, method = mycor, use = "complete.obs")
		txt <- format(c(r, 0.123456789), digits=digits)[1]
		txt <- paste(prefix, txt, sep="")
		if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
		text(0.5, 0.5, txt, cex = cex.cor^abs(r))
		}

smooth.type <- match.arg(smooth.fun)
switch(smooth.type, LM = (myfun = mylm), RLM = (myfun = myrlm), lmRob = (myfun = lmRob), SMOOTH = (myfun = panel.smooth))		
cor.type <- match.arg(cor.method)
switch(cor.type, spearman = (mycor = "spearman"), pearson = (mycor = "pearson"))		
			
pairs(data,
	panel=function(x,y)myfun(x, y, Pcol, Lcol,...),
	upper.panel=function(x, y)panel.cor(x, y, mycor = mycor, cex.cor = CexCor))
if(subTitle) title(sub = paste("Lower panel : Smoothing by", smooth.type, " / Upper panel : Correlation by", cor.type), ...)
}

