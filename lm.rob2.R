lm.rob2 <- function (formula, Plot=FALSE) {

		x <- x
		y <- y
		na.index <- which(is.na(x) | x==Inf | x==(-Inf) | is.na(y) | y==Inf | y==(-Inf))
		if(length(na.index)>0) {x <- x[-na.index]; y <- y[-na.index]
						cat(length(na.index), "values suppressed due to missingness\n")}

		epsilon = 1e-5
		c = 5

		lm1 <- lm(formula)					#; abline(lm1)
		s <- summary(lm1)$sigma
		s.init = Inf

		j = 0							# compteur itérations

		while((s.init - s) > 1e-10){

			s.init <- s

			# Pondérations
				r <- resid(lm1)
				m = median (r, na.rm = T)
				s = median (abs (r - m), na.rm = T)
				u = (r - m) / ((c * s) + epsilon)
				W = rep (0, length(x))
				i = abs(u) <= 1
				w <- rep(0, length(x))
				w[i] = ((1 - u^2)^2)[i]
			
			# régression pondérée
			lm1 <- lm(formula, weights = w)			#; abline(lm1)

			# nouveaux résidus
			s <- summary(lm1)$sigma
			j = j+1
			# print(paste("iter =", j))

			}

		if (Plot){
			plot(x,y)
			abline (lm1, col="red")
			}

		return(lm1)
}		

