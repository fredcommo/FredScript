ComBat3 <- function(eset, info, type='txt', write=F, covariates= NULL, par.prior=T, filter=F, prior.plots=F){

# info : contains experiments annotaions. 
#	the 1st column must contains expIds identical to eset colnames.
#	the 2nd column must be named 'Batch' and contain the batch annotations

# covariates : which covariates have to be taken into account.
# 	Key the specific col numbers if needed, except 1:2. ex : covariates = c(3, 5)
#	Key NULL to take into account only SampleIds and Batches (default)
#	Key 'all' to take into account all the available covariates

# Additional functions

	# filters data based on presence/absence call
	filter.absent <- function(x,pct){
		present <- T
		col <- length(x)/2
		pct.absent <- (sum(x[2*(1:col)]=="A") + sum(x[2*(1:col)]=="M"))/col
		if(pct.absent > pct){present <- F}
		present
	}

	# Next two functions make the design matrix (X) from the sample info file 
	build.design <- function(vec, des=NULL, start=2){
		tmp <- matrix(0,length(vec),nlevels(vec)-start+1)
		for (i in 1:ncol(tmp)){tmp[,i] <- vec==levels(vec)[i+start-1]}
		cbind(des,tmp)
	}

	design.mat <- function(saminfo){
		tmp <- which(colnames(saminfo) == 'Batch')
		tmp1 <- as.factor(saminfo[,tmp])
		cat("Found",nlevels(tmp1),'batches\n')
		design <- build.design(tmp1,start=1)
		ncov <- ncol(as.matrix(saminfo[,-c(1:2,tmp)]))
		cat("Found", ncov, 'suppl. covariate(s)\n')
		if(ncov>0){
			for (j in 1:ncov){
				tmp1 <- as.factor(as.matrix(saminfo[,-c(1:2,tmp)])[,j])
				design <- build.design(tmp1,des=design)
				}
			}
		design
	}

	# Makes a list with elements pointing to which array belongs to which batch
	list.batch <- function(saminfo){
		tmp1 <- as.factor(saminfo[,which(colnames(saminfo) == 'Batch')])
		batches <- NULL
		for (i in 1:nlevels(tmp1)){batches <- append(batches, list((1:length(tmp1))[tmp1==levels(tmp1)[i]]))}
		batches
	}

	# Trims the data of extra columns, note your array names cannot be named 'X' or start with 'X.'
	trim.dat <- function(eset){
		tmp <- strsplit(colnames(eset),'\\.')
		tr <- NULL
		for (i in 1:length(tmp)){tr <- c(tr,tmp[[i]][1]!='X')}
		tr
	}

	# Following four find empirical hyper-prior values
	aprior <- function(gamma.hat){m=mean(gamma.hat, na.rm = T); s2=var(gamma.hat, na.rm = T); (2*s2+m^2)/s2}
	bprior <- function(gamma.hat){m=mean(gamma.hat, na.rm = T); s2=var(gamma.hat, na.rm = T); (m*s2+m^3)/s2}
	postmean <- function(g.hat,g.bar,n,d.star,t2){(t2*n*g.hat+d.star*g.bar)/(t2*n+d.star)}
	postvar <- function(sum2,n,a,b){(.5*sum2+b)/(n/2+a-1)}


	# Pass in entire data set, the design matrix for the entire data, the batch means, the batch variances, priors (m, t2, a, b), columns of the data  matrix for the batch. Uses the EM to find the parametric batch adjustments

	it.sol  <- function(sdat,g.hat,d.hat,g.bar,t2,a,b,conv=.0001){
		n <- apply(!is.na(sdat), 1, sum, na.rm = T)
		g.old <- g.hat
		d.old <- d.hat
		change <- 1
		count <- 0
		while(change>conv){
			g.new <- postmean(g.hat,g.bar,n,d.old,t2)
			sum2 <- apply((sdat-g.new%*%t(rep(1,ncol(sdat))))^2, 1, sum, na.rm = T)
			d.new <- postvar(sum2,n,a,b)
			change <- max(abs(g.new-g.old)/g.old,abs(d.new-d.old)/d.old, na.rm = T)
			g.old <- g.new
			d.old <- d.new
			count <- count+1
			}
		#cat("This batch took", count, "iterations until convergence\n")
		adjust <- rbind(g.new, d.new)
		rownames(adjust) <- c("g.star","d.star")
		adjust
	}

	#likelihood function used below
	L <- function(x,g.hat,d.hat){prod(dnorm(x,g.hat,sqrt(d.hat)))}

	# Monte Carlo integration function to find the nonparametric adjustments
	int.eprior <- function(sdat,g.hat,d.hat){
		g.star <- d.star <- NULL
		r <- nrow(sdat)
		for(i in 1:r){
			g <- g.hat[-i]
			d <- d.hat[-i]		
			x <- sdat[i,!is.na(sdat[i,])]
			n <- length(x)
			j <- numeric(n)+1
			eset <- matrix(as.numeric(x),length(g),n,byrow=T)
			resid2 <- (eset-g)^2
			sum2 <- resid2%*%j
			LH <- 1/(2*pi*d)^(n/2)*exp(-sum2/(2*d))
			LH[LH=="NaN"]=0
			g.star <- c(g.star,sum(g*LH, na.rm = T)/sum(LH, na.rm = T))
			d.star <- c(d.star,sum(d*LH, na.rm = T)/sum(LH, na.rm = T))
			#if(i%%1000==0){cat(i,'\n')}
			}
		adjust <- rbind(g.star,d.star)
		rownames(adjust) <- c("g.star","d.star")
		adjust	
	} 

#fits the L/S model in the presence of missing data values

	Beta.NA <- function(y,X){
		des=X[!is.na(y),]
		y1=y[!is.na(y)]
		B <- solve(t(des)%*%des)%*%t(des)%*%y1
		B
	}
# End additional functions

# Start ComBat3

	tmp <- match(colnames(eset), info[,1])
	if(any(is.na(tmp))){stop('ERROR: Sample Information File and Data Array Names are not the same!')}
	tmp1 <- match(info[,1],colnames(eset))
	saminfo <- info[tmp1[!is.na(tmp1)],]		

	if(any(covariates != 'all')){saminfo <- saminfo[,c(1:2,covariates)]}
	if(is.null(covariates)){saminfo <- saminfo[,c(1:2)]}
	design <- design.mat(saminfo)	

	batches <- list.batch(saminfo)
	n.batch <- length(batches)
	n.batches <- sapply(batches, length)
	n.array <- sum(n.batches)
	

	## Check for missing values
	NAs = any(is.na(eset))
	if(NAs){cat(c('Found',sum(is.na(eset)),'Missing Values in eset\n'),sep=' ')}


	##Standardize Data across genes
	cat('Standardizing eset across genes\n')
	if (!NAs){B.hat <- solve(t(design)%*%design)%*%t(design)%*%t(as.matrix(eset))} else{B.hat=apply(eset,1,Beta.NA,design)} #Standarization Model
	grand.mean <- t(n.batches/n.array)%*%B.hat[1:n.batch,]
	if (!NAs){var.pooled <- ((eset-t(design%*%B.hat))^2)%*%rep(1/n.array,n.array)} else{var.pooled <- apply(eset-t(design%*%B.hat),1,var,na.rm=T)}

	stand.mean <- t(grand.mean)%*%t(rep(1,n.array))
	if(!is.null(design)){tmp <- design;tmp[,c(1:n.batch)] <- 0;stand.mean <- stand.mean+t(tmp%*%B.hat)}	
	s.data <- (eset-stand.mean)/(sqrt(var.pooled)%*%t(rep(1,n.array)))

	##Get regression batch effect parameters
	cat("Fitting L/S model and finding priors\n")
	batch.design <- design[,1:n.batch]
	if (!NAs){gamma.hat <- solve(t(batch.design)%*%batch.design)%*%t(batch.design)%*%t(as.matrix(s.data))}else{gamma.hat=apply(s.data,1,Beta.NA,batch.design)}
	delta.hat <- NULL
	for (i in batches){
		delta.hat <- rbind(delta.hat,apply(s.data[,i], 1, var, na.rm=T))
		}

	##Find Priors
	gamma.bar <- apply(gamma.hat, 1, mean, na.rm = T)
	t2 <- apply(gamma.hat, 1, var, na.rm = T)
	a.prior <- apply(delta.hat, 1, aprior)
	b.prior <- apply(delta.hat, 1, bprior)

	
	##Plot empirical and parametric priors

	if (prior.plots & par.prior){
		par(mfrow=c(2,2))
		tmp <- density(gamma.hat[1,])
		plot(tmp,  type='l', main="Density Plot")
		xx <- seq(min(tmp$x, na.rm = T), max(tmp$x, na.rm = T), length=100)
		lines(xx,dnorm(xx,gamma.bar[1],sqrt(t2[1])), col=2)
		qqnorm(gamma.hat[1,])	
		qqline(gamma.hat[1,], col=2)	
	
		tmp <- density(delta.hat[1,])
		invgam <- 1/rgamma(ncol(delta.hat),a.prior[1],b.prior[1])
		tmp1 <- density(invgam)
		plot(tmp,  typ='l', main="Density Plot", ylim=c(0, max(tmp$y, tmp1$y, na.rm = T)))
		lines(tmp1, col=2)
		qqplot(delta.hat[1,], invgam, xlab="Sample Quantiles", ylab='Theoretical Quantiles')	
		lines(c(0,max(invgam)),c(0,max(invgam)),col=2)	
		title('Q-Q Plot')
	}
	
	##Find EB batch adjustments

	gamma.star <- delta.star <- NULL
	{
	if(par.prior){
		cat("Finding parametric adjustments\n")
		for (i in 1:n.batch){
			temp <- it.sol(s.data[,batches[[i]]],gamma.hat[i,],delta.hat[i,],gamma.bar[i],t2[i],a.prior[i],b.prior[i])
			gamma.star <- rbind(gamma.star,temp[1,])
			delta.star <- rbind(delta.star,temp[2,])
			}
   		}
  	else{
		cat("Finding nonparametric adjustments\n")
		for (i in 1:n.batch){
			temp <- int.eprior(as.matrix(s.data[,batches[[i]]]),gamma.hat[i,],delta.hat[i,])
			gamma.star <- rbind(gamma.star,temp[1,])
			delta.star <- rbind(delta.star,temp[2,])
			}
		}
 	}

	### Normalize the Data ###
	cat("Adjusting the Data\n")

	bayesdata <- s.data
	j <- 1
	for (i in batches){
		bayesdata[,i] <- (bayesdata[,i]-t(batch.design[i,]%*%gamma.star))/(sqrt(delta.star[j,])%*%t(rep(1,n.batches[j])))
		j <- j+1
		}

	bayesdata <- (bayesdata*(sqrt(var.pooled)%*%t(rep(1,n.array))))+stand.mean

	{
	if(write){
		output_file <- paste('Adjusted',eset,'.xls',sep='_')
                 #print(geneinfo[1:2])
                 #print(bayesdata[1:2,1:4])
		 #cat(c(colnames(geneinfo),colnames(dat),'\n'),file=output_file,sep='\t')
		#suppressWarnings(write.table(cbind(geneinfo,formatC(as.matrix(bayesdata), format = "f")), file=output_file, sep="\t", quote=F,row.names=F,col.names=F,append=T))
    outdata <- cbind(ProbeID=geneinfo, bayesdata); write.table(outdata, file=output_file, sep="\t")
		cat("Adjusted data saved in file:",output_file,"\n")
		}
    else{return(bayesdata)} # else{return(cbind(geneinfo,bayesdata))}
	}

}
