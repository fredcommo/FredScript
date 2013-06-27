IC50.5P <- function(dose, resp, sd = NULL, Wcoef = 0.25, Plot = TRUE, unit="nM", Title="", Plot.sd = TRUE, AddIC = TRUE, Add = FALSE, output = TRUE, pCol = "grey50", lCol = "grey25", subTit = TRUE,...){

	# 5P logistic function
		L5P <- function(bottom, top, xmid, scal, s,  X){
			Y.fit <-bottom+(top-bottom)/(1+10^((xmid-X)*scal))^s
			return(Y.fit)
			}

	# Weighted SCE Function (Sum of Squared errors)
		sce.5P <- function(param, X, yobs, Weights, Wcoef) {
			bottom <- param[1]
			top <- param[2]
			xmid <- param[3]
			scal <- param[4]
			s <- param[5]
			ytheo <- L5P(bottom, top, xmid, scal, s, X)
			residus <- yobs - ytheo
			Weights <- (1/(residus^2))^(Wcoef)
			return(sum(Weights*(yobs - ytheo)^2))
			}

	# Get weights (not used)
		sce.5P.diag <- function(yobs, ytheo, w) {
			sq.res <- (yobs - ytheo)^2
			weights <- 1/sq.res^w
			return(weights)
			}

	# initialization of parameters
		# add.first = FALSE		
		#if(y[1]<0.9){
			x.marge <- (max(dose)-min(dose))/(length(dose)-1)
			x <- c(min(dose)-x.marge, dose)
			y <- c(1, resp)
			sd <- c(0, sd)
			add.first = TRUE
		#}

		bottom.ini <- miny <- min(y, na.rm = T)
		top.ini = maxy = max(y, na.rm = T)
		minx = min(x, na.rm = T)
		maxx = max(x, na.rm = T)
		xmid.ini = (maxx + minx)/2
		z <- (y - bottom.ini)/(top.ini - bottom.ini)
		z[z==0] <- 0.01; z[z==1] <- 0.99
		scal.ini = coef(lm(x ~ log(z/(1-z))))[2]
		scal.ini <- as.numeric(scal.ini); scal.ini
		s.ini = 1
		weights <- rep(1, length(y))		

	# Optimisation step using nlm()
		init <- c(bottom.ini, top.ini, xmid.ini, scal.ini, s.ini)
		best <- nlm(f = sce.5P, p = init, X = x, yobs = y, Weights = weights, Wcoef = Wcoef)
		
	# Get best parameters
		best.bottom <- best$estimate[1]
		best.top <- best$estimate[2]
		best.xmid<-best$estimate[3]
		best.scal <- best$estimate[4]
		best.s <- best$estimate[5]

	# Goodness of fit
		yfit <- L5P(best.bottom, best.top, best.xmid, best.scal, best.s, dose)
		#weights <- sce.5P.diag(resp, yfit, Wcoef)
		lm.test <- lm(yfit ~ resp)	#, weights = weights)
		Goodness <- summary(lm.test)
		
	# Estimate critical points
		newX <- seq(minx - x.marge, maxx + x.marge, length = 100)						
		Y.best <- L5P(best.bottom, best.top, best.xmid, best.scal, best.s, newX)

			# Inflexion point using the second derivative
				Xflex = best.xmid + (1/best.scal)*log10(best.s)
				Yflex = best.bottom + (best.top - best.bottom)*(best.s/(best.s+1))^best.s

			# IC50 x.coordinate
				Y50 = 0.5
				X50 = best.xmid - 1/best.scal*log10(((best.top - best.bottom)/(Y50 - best.bottom))^(1/best.s)-1)

			# Slope at inflexion point using the first derivative
				B = best.bottom + (best.top - best.bottom)*log(10)*(best.scal)*(best.s/(best.s+1))^(best.s+1)
				A = Yflex  - B*(Xflex)

			# Compute simulations to estimate the IC50 conf. interval
				ytheo <- L5P(best.bottom, best.top, best.xmid, best.scal, best.s,  x)
				Q <- sum((y - ytheo)^2)
				S2 <- Q/(length(x)-2)
				s <- sqrt(S2)

				if(is.na(X50)) min.D50 <- max.D50 <- NA
				else{
					MC.X50 <- rep(0, 10000)
					for (i in 1:10000){
						yi <- Y50 + rnorm(1, 0, s)
						MC.X50[i] <- best.xmid - 1/best.scal*log10(((yi-best.bottom)/(best.top-best.bottom))^(-1/best.s) - 1)
						}
					min.X50 <- quantile(MC.X50, probs=0.025, na.rm=T)
					min.D50 <- signif(10^min.X50, 2)
					max.X50 <- quantile(MC.X50, probs=0.975, na.rm=T)
					max.D50 <- signif(10^max.X50, 2)
					}

	# Graphics
	if(Plot){
		if(!Add) plot(y[-1] ~ x[-1], pch = 20, cex = 2.5, col = pCol,
							xlim = range(minx, maxx), ylim = range(0, 1),...)
		else points(y[-1] ~ x[-1], pch = 20, cex = 2.5, col = pCol, xlim = range(minx, maxx), ylim = range(0, max(maxy, Y.best)))

		if(AddIC){
			legend1 <- paste("IC50 :", round(10^X50, 2),unit)
			legend2 <- paste0("[", min.D50, " , ", max.D50, "]")
			legend('bottomleft', legend = c(legend1, legend2), cex = 1.2, text.col = 'steelblue', bty = 'n')
			}

		if(Plot.sd){
			pas = x.marge/4
			for(s in 2:length(sd)){
				segments(x0 = x[s], x1 = x[s], y0 = y[s]-sd[s], y1 = y[s]+sd[s], lty = 3)
				segments(x0 = x[s]-pas, x1 = x[s]+pas, y0 = y[s]-sd[s], y1 = y[s]-sd[s], lty = 3)
				segments(x0 = x[s]-pas, x1 = x[s]+pas, y0 = y[s]+sd[s], y1 = y[s]+sd[s], lty = 3)
				}
			}

		lines(Y.best ~ newX, col = lCol,...)
		Sub = NULL
		if(subTit){
			Sub = "Weighted 5P logistic regr. (F. Commo)"
			if(Wcoef==0) Sub = "Non weighted 5P logistic regr. (F. Commo)"
			}
		title (main = Title, sub = Sub)
		}

		# Output
		if(output){
		return(list(init = init, optim = cbind.data.frame(bottom = best.bottom, top = best.top, xmid = best.xmid, scal = best.scal, s = best.s),
					Xflex = Xflex, Yflex = Yflex, X50 = X50, Slope = B, slopeIntercept = A,
					D50 = 10^X50, Dmin = min.D50, Dmax = max.D50,
					xFit = newX, yFit = Y.best, Goodness = Goodness))
		}
}
