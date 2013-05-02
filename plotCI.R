plotCI<-function (y, datax = FALSE, p=0.05,...) 
{
    
	my<-mean(y)
	y <- quantile(y[!is.na(y)], c(0.25, 0.75))
	x <- qnorm(c(0.25, 0.75))
	x1 <- quantile(y[!is.na(y)],p)
	x2<- quantile(y[!is.na(y)],1-p)

    	if (datax) {
    		slope <- diff(x)/diff(y)
    		}

    	else {
        	slope <- diff(y)/diff(x)
     		}
    	abline(my,slope)
	abline(x1, slope,lty=3,col="red", ...)
	abline(x2, slope,lty=3,col="green", ...)
}
