plot.freq.grad <- function (xobs , yobs , ncla = 10, showcla=T, scale.x=F, scale.y=F, Plot=T, title=NULL, xlab=NULL, ylab="Proba", na.rm=TRUE) 

# xobs 	:	ordered x-values (numeric)
# yobs 	:	a success/failure vector, ordered according to x
# ncla	:	a number of intervals
# showcla	:	do the intervals have to be drawn?
# scale.x	:	If TRUE, x-values will be normalized
# scale.y	:	If TRUE, y-range will be forced to [0,1]

{
    if(any(is.na(xobs), is.na(yobs))){
	NAs <- which(is.na(xobs) | is.na(yobs))
	xobs <- xobs[-NAs]
	yobs <- yobs[-NAs]
	}

    if (is.numeric(xobs) == F) return("xobs:numeric expected")
    if (scale.x) xobs <- scale(xobs)

    yobs <- yobs[order(xobs)]
    xobs <- sort(xobs)
    xobs <- xobs + rnorm(length(xobs), 0, 1e-4)	

    q0 <- quantile(xobs, probs = seq(0, 1, 1/ncla), na.rm = na.rm)
 
    # c.cla <- quantile(xobs, probs = seq(0 + 1/2/ncla, 1 - 1/2/ncla, 1/ncla), na.rm = na.rm)

    c.cla <- rep(0, length(q0)-1)
    for (i in 1:(length(q0)-1)) c.cla[i] <- (q0[i]+q0[i+1])/2

    xmin <- min(xobs, na.rm=T) ; xmax <- max(xobs,na.rm=T)
    q1 <- cut(xobs, q0, include.lowest = T)

    t0 <- table(q1, yobs)
    # t0[, 1] <- (t0[, 1] + t0[, 2])
    freq <- t0[, 2]/rowSums(t0)

    basbar <- rep(0, ncla) ; haubar <- rep(0, ncla)
    
	for (i in 1:ncla) {
        succes <- t0[i, 2] ;  essai <- sum(t0[i, ])
        if (essai >= 8) {
            a0 <- prop.test(succes, essai)$conf.int
            basbar[i] <- a0[1] ; haubar[i] <- a0[2]
            if (a0[1] > freq[i]) basbar[i] <- NA
            if (a0[2] < freq[i]) haubar[i] <- NA
        	} 

		else {
            basbar[i] <- NA ; haubar[i] <- NA
        	}
    	}

    ymin <- min(basbar, na.rm = T)
    ymin <- min(ymin,min(freq))
    ymax <- max(haubar, na.rm = T)
    ymax <- max(ymax,max(freq))    

    if (scale.y) { ymin=0; ymax=1}

    if(Plot){
    	plot(xobs, yobs, ylim = range(0, 1), xlim = c(xmin, xmax), xlab=xlab, ylab = ylab, type = "n", main=title)
    	points(c.cla, freq, pch = 16,col="red")
	

    	for (i in 1:ncla) {
      	if (!is.na(basbar[i])) {
            	size.bar <- par("cxy")[1]/4
            	segments(c.cla[i], basbar[i] , c.cla[i], haubar[i],col="blue" )
            	segments(c.cla[i] - size.bar, haubar[i] , c.cla[i] + size.bar, haubar[i],col="blue" )
            	segments(c.cla[i] - size.bar, basbar[i] , c.cla[i] + size.bar, basbar[i],col="blue" )
        		}
        	if ((i > 1) & showcla) abline(v = q0[i], lty = 2)
    		}
	}

	return(list(centres=c.cla, freq=freq, freq.tab=t0))
}
 

