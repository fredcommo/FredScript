plot.logistic <- function(x, y, Xlab = "", Ylab = "Probability of class", Plot.Freq = TRUE, Ncla = 10, Title = NULL, add.p = TRUE, Col = c("royalblue3","indianred2")){

# x : values
# y : vector of classes
# Plot.Freq : add frequences by interval
# Ncla : define the number of intervals
# add.p : add the p-value for the logistic model

	source("D:\\Stats\\Doc R\\Scripts R\\plot.freq.grad.R")

	y <- as.factor(y[order(x)])
	x <- sort(x)
	
	model <- glm(y~x, family=binomial)
	
	if (Plot.Freq){
		plot.freq.grad(xobs = x, yobs = y, showcla = F, ncla = Ncla, scale.x = F, xlab = Xlab, ylab = Ylab, title = Title)
		points(c(0,1)[y]~x, pch = c(17,19)[y], col=Col[y], cex=1.5)
		}

	else plot(c(0,1)[y]~x, pch = c(17,19)[y], col=Col[y], cex=1.5, xlab = Xlab, ylab = Ylab)

	newx <- seq(min(x), max(x), len=100)
	lines(predict(model, newdata=list(x=newx), type="response") ~ newx, col = "cyan", lwd = 3)
	axis(side = 4, at = c(0,1), labels = levels(y), font = 4, cex=1.5)

	if(add.p){
		slope <- model$coefficients[2]
		pval <- summary(model)$coefficients[2,4]
		xpos <- min(x)*0.9 + 0.1*max(x)
		ypos <- 0.8
		if(slope<0) ypos <- 0.2 
		text(x = xpos, y = ypos, labels=paste("p-value", "\n",signif(pval, 3)), cex = 1.2, font = 4, col = "darkblue")
		}

	print(summary(model))
}


