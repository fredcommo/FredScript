
model2P <- function(x, xmid, scal){
	y.fit <- 1/(1+exp((xmid-x)/scal))
	return(y.fit)
	}
		

C <- function(x, xmid, scal){
	Cx <- 1/scal*(x - xmid)
	return(Cx)
	}

P <- function(x, xmid, scal){
	p <- 1/(1+exp(-C(x, xmid, scal)))
	return(p)
	}

dLL <- function(par){
	xmid <- par[1]
	scal <- par[2]	
	dll <- (-2)*sum(y*log(P(x, xmid, scal))+(1-y)*log(1-P(x, xmid, scal)))
	return(dll)
	}

# glm.res <- data.frame(row.names=rownames(comb), Index=seq(1,nrow(comb)), Gene=NA, odds=0, pval=0, slope=0)


par(mfrow=c(2,4))
index = as.vector(which(Symb=="CDH11"))

# glm.res <- data.frame(row.names=rownames(comb), Index=seq(1,nrow(comb)), Gene=NA, odds=0, pval=0, slope=0)

for (g in index){
	gname <- gettext(Symb[g])
	probe <- names(Symb[g])

	x <- as.numeric(eset[g, select])
	y <- ifelse(cad.grp=="hi.cdh2", 1, 0)

	xmid.ini = (max(x) + min(x))/2
	scal.ini = 1/coef(lm1 <- lm(x~y))[2]

	opt <- optim(c(xmid = xmid.ini, scal = 1/scal.ini), dLL, hessian=F)

	best.xmid <- opt$par[1]
	best.scal <- opt$par[2]

	newx <- seq(min(x), max(x), len=100)
	post <- P(newx, best.xmid, best.scal)					# model2P(newx, best.xmid, best.scal)

	plot(y ~ x, pch=c(17,19)[cad.grp], col=c("deepskyblue","violetred4")[cad.grp], cex=1)
	lines(post ~ newx, lwd=2, col="orangered")
	axis(side=4, at=c(0,1), labels=c("Hi.CAD11","Hi.NCAD"))
	abline(v=median(x[y==0]), lty=3, col="deepskyblue")
	abline(v=median(x[y==1]), lty=3, col="violetred4")
	abline(v=best.xmid, lty=3, col="red")
	abline(h=0.5, lty=3, col="red")

	summary(glm1 <- glm(y~x, family=binomial))
	odds <- summary(glm1)$coefficients[2,1]
	pval <- summary(glm1)$coefficients[2,4]
	pente = 1/best.scal; pente

	BoxPoints(cad.grp, x, disp = 0.02)	#, col=c("deepskyblue","violetred4"))
	abline(h=best.xmid, lty=3, col="red")
	title(main = paste(gname, ";", probe))

	# glm.res[g, ]<- c(g, gname, signif(odds,4), signif(pval,4), signif(pente,4))
	# cat(g, "of", nrow(comb),"\n")
	}

	# pred <- predict(glm1, newdata=list(x=newx), type="response")
	# lines(pred ~ newx, col="forestgreen")
	  title(main = paste(gname, ";", probe))
	# text(x = max(x)*0.75, y = 0.5, labels = paste("Odds = ", signif(odds, 3), ",", "p =", signif(pval,3)),
	#		font = 3, col="violetred4")


