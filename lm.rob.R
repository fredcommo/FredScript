lm.rob <- function (x, y, Eps = 1e-5, C = 5, Plot=FALSE) {

		na.index <- which(is.na(x) | x==Inf | x==(-Inf) | is.na(y) | y==Inf | y==(-Inf))
		if(length(na.index)>0) {x <- x[-na.index]; y <- y[-na.index]
						cat(length(na.index), "values suppressed due to missingness\n")}

		epsilon = Eps
		c = C

		lm.test <- lm(y ~ x)					#; abline(lm1)
		s <- summary(lm.test)$sigma
		s.init = Inf

		j = 0							# compteur iterations

		while((s.init - s) > 1e-10){

			s.init <- s

			# Ponderations
				r <- resid(lm.test)
				m = median (r, na.rm = T)
				s = median (abs (r - m), na.rm = T)
				u = (r - m) / ((c * s) + epsilon)
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

		if (Plot){
			plot(x,y)
			abline (lm.test, col="red")
			}

		return(list(model = lm.test, x = x, y = y, W = w))
}		

