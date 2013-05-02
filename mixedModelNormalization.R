require(nlme)

mixedModelNormalization <- function(log2R,log2G,array=0,geneID=0, row=0, column=0)
{
if(length(array) <= 1)
	array <- rep(1,length(log2R))
	
#Count the number of arrays
u.array<-unique(array)
n.array<-length(u.array)
aics <- matrix(999999999,nrow=n.array,ncol=7)

par(mfrow=c(1,1))
for (i in sort(u.array))
{	
	############################################
	## Select the data from the current array ##
	############################################
	
	log2Ri <- as.vector(log2R[array==u.array[i]])
	log2Gi <- as.vector(log2G[array==u.array[i]])
	if(length(geneID) > 1)
		geneIDi  <- as.vector(geneID[array==u.array[i]])

	###################
	## Sort the data ##
	###################
	Ax<-(log2Ri+log2Gi)/2
	Mx<-log2Ri-log2Gi

	A<-sort(Ax)
	ord <- order(Ax)
	M<-Mx[ord]

	ny<-length(A)
	nx1<-20
		
	if(length(geneID) > 1)
	{
		geneIDi <- geneIDi[ord]
		all <-rep(1,length(geneIDi))
	}
	else
	{
		all <- rep(1,length(A))
	}
	
	if(length(row) > 1)
	{
	
	################################
	## Find the resp. pin numbers ##
	################################
	
	printTip <-rep(1,length(A))
	if(length(column) > 1)
	{
		rowi <- row[array==u.array[i]]
		columni<-column[array==u.array[i]]
		maxcol <- max(columni)
		for (j in 1:length(A))
			printTip[j]<-((rowi[j]-1)%%maxcol)*maxcol + columni[j]
	}
	else
	{
		printTip <- row
	}
		

	#############################
	## Convert to RI-plot data ##
	#############################

	printTip <- printTip[ord]
	nrOfPins <- length(unique(printTip))

	#############################################
	## Make "nrOfPins" groups (for pin effect) ##
	#############################################

	pins <- matrix(0,ny,nrOfPins-1)
	for (j in 1:ny)
		if(printTip[j] < nrOfPins)
			pins[j,printTip[j]] = 1
}
	#######################################
	## Chop x in pieces for the z-matrix ##
	#######################################


	knots <- quantile(A,probs=seq(from=0.000000000000001,to=1-0.000000000000001,length=nx1))
	zmat1<-matrix(0,ny,nx1)	

	for(j in 1:nx1)
		for(k in 1:ny)
			if(A[k]>knots[j])
				zmat1[k,j]<-A[k]-knots[j]

	##################################
	## Fit several different models ##
	##################################

	lmm.global <- lme(M ~ 1,random=~1|all)
	lmm.linear <- lme(M ~ A,random=~1|all)
	lmm.nonlin <- lme(M ~ A, random = list(all = pdIdent(~zmat1)))

	if(length(geneID) > 1)
		lmm.randint<-lme(M ~ A,random=list(all=pdIdent(~zmat1-1),geneIDi = pdSymm(~1)))

	if(length(row) > 1)
	{
		lmm.fixed <- lme(M ~ A  + pins, random = list(all=pdIdent(~zmat1)))
		lmm.pin <-lme(M ~ A,random=list(printTip=pdIdent(~zmat1 -1)))
		if(length(geneID) >1)
			lmm.randint.pin <-lme(M ~ A,random=list(printTip=pdIdent(~zmat1 -1),geneIDi = pdSymm(~1)))
	}
	
	
	######################################################
	## Compute the AICs and determine which is the best ##
	######################################################
	
	aics[i,1] <- AIC(lmm.global)
	aics[i,2] <- AIC(lmm.linear)
	aics[i,3] <- AIC(lmm.nonlin)
	if(length(geneID) > 1)
		aics[i,4] <- AIC(lmm.randint)
	if(length(row) > 1)
	{
		aics[i,5] <- AIC(lmm.fixed)
		aics[i,6] <- AIC(lmm.pin)
		if(length(geneID) > 1)
			aics[i,7] <- AIC(lmm.randint.pin)
	}
	print(aics[i,])	
	o <- order(aics[i,])

	if(o[1] == 1)
	{
		y.pred <- as.vector(predict(lmm.global))
				
	}
	else if(o[1] == 2)
	{
		y.pred <- as.vector(predict(lmm.linear))

	}
	else if(o[1] == 3)
	{
		b.vec1<-lmm.nonlin$coeff$random$all
		zb1<-zmat1%*%t(b.vec1)
		beta.vec1<-lmm.randint$coeff$fixed
		xmat1<-cbind(rep(1,length(A)),A)
		xbeta1<-xmat1%*%t(beta.vec1)
		y.pred <- xbeta1+zb1
	}
	else if(o[1] == 4)
	{
		b.vec2<-lmm.randint$coeff$random$all
		zb2<-zmat1%*%t(b.vec2)
		beta.vec2<-t(lmm.randint$coeff$fixed)
		xmat1<-cbind(rep(1,length(A)),A)
		xbeta2<-xmat1%*%t(beta.vec2)
		y.pred <- xbeta2+zb2
	}
	else if(o[1] == 5)
	{	
		b.vec3 <- lmm.fixed$coeff$rand$all
		zb3<-zmat1%*%b.vec3[-1]
		beta.vec3<-lmm.fixed$coeff$fixed[1:2]
		xmat1<-cbind(rep(1,length(A)),A)
		xbeta3<-xmat1%*%beta.vec3
		pins3 <- lmm.fixed$coeff$fixed[3:(nrOfPins+1)]
		prod <- pins %*% pins3
	
		y.pred <- xbeta3 + zb3 + prod
	}
	else if(o[1] == 6)
	{
		b.vec4<-lmm.pin$coeff$rand$printTip
		zb4<-zmat1%*%t(b.vec4)
		beta.vec4<-t(lmm.pin$coeff$fixed)
		xmat1<-cbind(rep(1,length(A)),A)
		xbeta4<-xmat1%*%t(beta.vec4)
		y.pred <- rep(0,length(A))
		for(j in 1:length(A))
			y.pred[j] <- xbeta4[j]+zb4[j,printTip[j]]
		
	}	
	else if(o[1] == 7)
	{
		b.vec3<-lmm.randint.pin$coeff$rand$printTip
		zb3<-zmat1%*%t(b.vec3)
		beta.vec3<-t(lmm.randint.pin$coeff$fixed)
		xmat1<-cbind(rep(1,length(A)),A)
		xbeta3<-xmat1%*%t(beta.vec3)
		y.pred <- rep(0,length(A))
		for(j in 1:length(A))
			y.pred[j]<-xbeta3[j]+zb3[j,printTip[j]] 



	}
	z <- M - y.pred
	windows()
	plot(A, M,xlab= "A", ylab="M",cex=1.5,pch=".")
	title(paste("Array", as.character(i), "before normalization"))
	
	if(o[1] < 5)
	{
		lines(A,y.pred)
	}
	else
	{
		for (pin in 1:nrOfPins)
		{	print(printTip)
			Atemp <- A[printTip == pin]
			Asorted <- sort(Atemp)
			Mtemp <- y.pred[printTip == pin]
			ord <- order(Atemp)
			Msorted <- Mtemp[ord]
			print(length(Asorted))
			print(length(Msorted))			
			lines(Atemp,Mtemp,lwd=2,col=pin)
		}
	}

	windows()
	plot(A,z,pch=".",xlab= "A", ylab="M",cex=1.5)			
	title(paste("Array", as.character(i), "after normalization"))
	
	log3na <- rep(0,length(A))
	log5na <- rep(0,length(A))

	log3na[ord] <- (A + 1/2*z)
	log5na[ord] <- (A - 1/2*z)

}

result <- matrix(0,nrow=length(log3na),ncol=2)
result[,1] <- log3na 
result[,2] <- log5na

invisible(result)

}
