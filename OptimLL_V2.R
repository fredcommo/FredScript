
model2P <- function(x, A, xmid, scal){
	y.fit <- (A)/(1+exp((xmid-x)/scal))
	return(y.fit)
	}
		

C <- function(x, xmid, scal){
	Cx <- 1/(scal*(x - xmid))
	return(Cx)
	}

P <- function(x, A, xmid, scal){
	p <- A/(1+exp(-C(x, xmid, scal)))
	return(p)
	}

dLL <- function(par){
	A <- par[1]
	xmid <- par[2]
	scal <- par[3]
	dll <- (-2)*sum(y*log(P(x, A, xmid, scal))+(1-y)*log(1-P(x, A, xmid, scal)))
	return(dll)
	}

par(mfrow=c(1,1))
g = 8473

# glm.res <- data.frame(row.names=rownames(comb), Index=seq(1,nrow(comb)), Gene=NA, odds=0, pval=0, slope=0)

# for (g in 1:nrow(comb)){
	gname <- gettext(Symb[filtr[g]])
	probe <- names(Symb[filtr[g]])

	x <- as.numeric(eset.filtr[g,])
	y <- 	ifelse(cad==1, 0, 1) #ifelse(info$CAD=="hi.cdh2", 1, 0)

	xmid.ini = (max(x) + min(x))/2
	scal.ini = 1/coef(lm1 <- lm(x~y))[2]

	opt <- optim(c(A=1, xmid = xmid.ini, scal = scal.ini), dLL, hessian=F)

	best.A <- opt$par[1]
	best.xmid <- opt$par[2]
	best.scal <- opt$par[3]

	newx <- seq(min(x), max(x), len=100)
	post <- model2P(newx, best.A, best.xmid, best.scal)					# model2P(newx, best.xmid, best.scal)

	plot(y ~ x, pch=c(17,19)[cad], col=c("deepskyblue","violetred4")[cad], cex=1.5)
	lines(post ~ newx, lwd=2, col="orangered")
	axis(side=4, at=c(0,1), labels=c("Hi.CAD11","Hi.NCAD"))

	summary(glm1 <- glm(y~x, family=binomial))
	odds <- summary(glm1)$coefficients[2,1]
	pval <- summary(glm1)$coefficients[2,4]
	pente = 1/best.scal; pente

	# glm.res[g, ]<- c(g, gname, signif(odds,4), signif(pval,4), signif(pente,4))
	# cat(g, "of", nrow(comb),"\n")
# }

	pred <- predict(glm1, newdata=list(x=newx), type="response")
	lines(pred ~ newx, col="forestgreen")
	title(main = paste(gname, ";", probe))
	text(x = max(x)*0.75, y = 0.5, labels = paste("Odds = ", signif(odds, 3), ",", "p =", signif(pval,3)),
			font = 3, col="violetred4")

	BoxPoints(cad, x)	#, col=c("deepskyblue","violetred4"))
	title(main = paste(gname, ";", probe))


