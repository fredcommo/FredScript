
Q <- qchisq(p = seq(0.05, 0.95, by = 0.1), df = 2)
sigma <- matrix(c(1, 0, 0, 1), 2, 2)

#S<-var(X)
#M<-mean(X)
x <- seq(min(X[,1]), abs(min(X[,1])), length = 100)
y <- seq(min(X[,2]), abs(min(X[,2])), length = 100)
sigmainv <- solve(S)
a <- sigmainv[1, 1]
b <- sigmainv[2, 2]
c <- sigmainv[1, 2]
z <- outer(x, y, function(x, y) (a * x^2 + b * y^2 + 2 * c * x * y))

plot(X, pch=19, cex=0.2, col="lightblue3")
contour(x, y, z, col = "blue4", levels = Q, 
	labels = seq(from = 0.05, to = 0.95, by = 0.1), add = T)
