Resamp.T <- function(X, B = 1000, trimmed = 0){

	na.replace <- function(x, trimmed){
		qx <- quantile(x, probs=c(trimmed, 1-trimmed), na.rm = T)
		x[which(x<qx[1] | x>qx[2])] <- NA
		return(x)
		}

	N <- ncol(X)
	M.ini <- 1/N*rowSums(X, na.rm = T)
	t <- xbar <- c()
	for (i in 1:B){
		newX <- X[,sample(1:N, replace = T)]
		if(trimmed>0){
			newX <- apply(newX, 1, na.replace, trimmed = trimmed)
			newX <- t(newX)
			}
		n <- ncol(newX)
		M <- 1/n*rowSums(newX, na.rm = T)
		xbar <- cbind(xbar, M)
		if(i%%10==0) cat("B:", i, "\t")
		if(i%%200==0) cat("\n")
		}
		cat("\n")
		Mboot <- 1/B*rowSums(xbar, na.rm = T)
		Biais <- Mboot - M.ini
		S2boot <- 1/(B-1)*rowSums((xbar-Mboot)^2, na.rm = T)
		SEboot <- sqrt(S2boot)
		t <- Mboot/SEboot

		return(list(Bbar = Mboot, Biais = Biais, SEboot = SEboot, Tvalue = t))
}


# X <- cap.comb[1:100,]
# test.100 <- Resamp.T(X, B=100)
# test.1000 <- Resamp.T(X, B=1000)
# biais.bar.100 <- mean(test.100$Biais); biais.sd.100 <- sd(test.100$Biais)
# biais.bar.1000 <- mean(test.1000$Biais); biais.sd.1000 <- sd(test.1000$Biais)


# par(mfrow = c(1,2))
# plot(density(test.100$Biais), main = "100 resamplings")
# legend("topleft", legend = c(paste("mean =", signif(biais.bar.100,3)), paste("sd =", signif(biais.sd.100,3))), bty = "n")
# plot(density(test.1000$Biais), main = "1000 resamplings")
# legend("topleft", legend = c(paste("mean =", signif(biais.bar.1000,3)), paste("sd =", signif(biais.sd.1000,3))), bty = "n")
# par(op)
