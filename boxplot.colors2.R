
boxplot.colors2<-function(X, Class=NULL, g=1, a=0.5, Start=0.6, End=0.16, title=NULL) {

# X : values in any format, or a vector of values if Class is defined (factor)
# g, a : color (rainbow) parameters in [0,1]

	if (is.null(Class)){	

		if(!is.data.frame(X)) X<-as.data.frame(X)

		X.med <- apply(X,2,median)
		col.palette = rainbow(length(X.med), start = Start, end = End, gamma=g, alpha=a)

		box.col= as.numeric(rank(X.med))
		boxplot(X, col=col.palette[box.col], main=title)
	}

	else {
		X <- as.numeric(X)
		col.palette = rainbow(nlevels(Class), start = Start, end = End, gamma=g, alpha=a)
		
		X.med <- rep(0, nlevels(Class))
		for (i in 1:nlevels(Class)) X.med[i] <- median(X[which(Class==levels(Class)[i])], na.rm=TRUE)

		box.col= as.numeric(rank(X.med))
		boxplot(X~Class, col=col.palette[box.col], main=title)
	}

}
