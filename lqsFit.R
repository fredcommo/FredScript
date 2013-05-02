lqsFit <- function(x, y){

# Type II linear regression

	# Calculate sums, SSq, ...

		n <- length(x)
		Sx <- sum(x)
		Sy <- sum(y)

		xbar <- mean(x, na.rm=T)
		ybar <- mean(y, na.rm=T)

		U <- x - xbar
		V <- y - ybar

		Suv <- cov(x,y)
		Su2 <- var(x)
		Sv2 <- var(y)

		sigx <- sqrt(Su2/(n-1))
		sigy <- sqrt(Sv2/(n-1))

	# Calculate parameters

		m = (Sv2 - Su2 + sqrt(((Sv2 - Su2)^2) + (4 * Suv^2)))/(2 * Suv)
		b = ybar - m * xbar
		r = Suv / sqrt(Su2 * Sv2)

		sm = (m/r) * sqrt((1 - r^2)/n);
		sb1 = (sigy - sigx * m)^2;
		sb2 = (2 * sigx * sigy) + ((xbar^2 * m * (1 + r))/r^2);
		sb = sqrt((sb1 + ((1 - r) * m * sb2))/n);

		return(list(Slope = m, y.intercept = b, r = r, sd.slope = sm, sd.y.intercept = sb))
}

