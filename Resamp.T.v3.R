Resamp.T.v3 <- function(X, class, B = 1000, trimmed = 0){

# compute bootstrap t-values

	na.replace <- function(x, trimmed){
		qx <- quantile(x, probs=c(trimmed, 1-trimmed), na.rm = T)
		x[which(x<qx[1] | x>qx[2])] <- NA
		return(x)
		}

	stat.test <- function(x, grp){t.test(x~grp)$statistic}

	class <- as.factor(class)
	classlab <- levels(class)

	if(any(is.na(class))){
		NAs <- which(is.na(class))
		class <- class[-NAs]
		X <- X[,-NAs]
		cat(length(NAs), "case(s) suppressed for missingness\n")
		}

	split.X <- split(as.data.frame(t(X)), class)
	X1 <- t(split.X[[1]])
	X2 <- t(split.X[[2]])
	n1 <- ncol(X1)
	n2 <- ncol(X2)

	# M.ini <- lapply(split.X, mean, na.rm = T)
	t <- boot1 <- boot2 <- S1 <- S2 <- t.boot <- c()

	for (i in 1:B){
		newX1 <- X1[,sample(1:n1, replace = T)]
		newX2 <- X2[,sample(1:n2, replace = T)]
		if(trimmed>0){
			newX <- apply(newX, 1, na.replace, trimmed = trimmed)
			newX <- t(newX)
			}
		M1 <- 1/n1*rowSums(newX1, na.rm = T)
		M2 <- 1/n2*rowSums(newX2, na.rm = T)
		S1 <- cbind(S1, apply(newX1, 1, var, na.rm = T))
		S2 <- cbind(S2, apply(newX2, 1, var, na.rm = T))
		boot1 <- cbind(boot1, M1)
		boot2 <- cbind(boot2, M2)
		Grp <- c(rep(1, n1), rep(2, n2))
		t.boot <- cbind(t.boot, apply(cbind(newX1, newX2), 1, stat.test, grp = Grp))

		# if(i%%10==0) 
			cat("B:", i, "\n")
		# if(i%%200==0) cat("\n")
		}
		cat("\n")

		Mboot1 <- 1/B*rowSums(boot1, na.rm = T)
		Mboot2 <- 1/B*rowSums(boot2, na.rm = T)
		# Biais <- Mboot - M.ini
		S2boot1 <- 1/B*rowSums(S1, na.rm = T)
		S2boot2 <- 1/B*rowSums(S2, na.rm = T)
		# SEboot <- sqrt(S2boot)
		t <- apply(abs(t.boot), 1, quantile, probs = 0.5)
		t <- ifelse(Mboot2 - Mboot1<0, -t, t)
		p <- pt(abs(t), df = (n1+n2-2), lower.tail = F)*2
		adj.p <- p.adjust(p, method = "BH")
		return(cbind.data.frame(M.boot1 = Mboot1, Var.boot1 = S2boot1, M.boot2 = Mboot2, Var.boot2 = S2boot2, Tvalue = t, p.value = p, adj.p.value = adj.p))
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
