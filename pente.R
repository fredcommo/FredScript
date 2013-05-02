
pente<-function(x, n=512){
	d<-density(as.numeric(x),na.rm=T,n=n)
	n<-length(d$x)
	slope<-rep(0,n-1)
	
	for (i in 6:n)
		slope[i]<-(d$y[i]-d$y[i-5])/(d$x[i]-d$x[i-5])
	
	slope<-zapsmall(slope,3) 	# A revoir !!!!!
	sign<-rep(0,length(slope)-1)

	for (i in 1:length(sign))
		sign[i]<-slope[i]*slope[i+1]
	
	list(x=d$x, y=d$y, slope=slope, sign=sign)
}
	

	g=7967

	p<-pente(data[g,])	

# modifier !! voir Rohr.data, ex: gène 5565 , 3414, un mode non détecté.
	sign<-rep(0,length(p$slope)-1)
	for (i in 2:length(p$slope))
		sign[i]<-p$slope[i]*p$slope[i+1]
	
	par(mfrow=c(1,2))
	plot(density(as.numeric(data[g,],na.rm=T)))
	plot(p$slope~p$x,type="l")

# modifier pour ne sélectionner que les pentes descendantes !!!
	abline(v=p$x[which(sign<0)])
	abline(h=0)



	all2<-all[-7138,]

	nb.mode<-rep(0,dim(all2)[1])
	for (i in 1:dim(all2)[1])
	{	p<-pente(all2[i,])
		nb.mode[i]<-length(which(p$sign<0))
	}


Duke.data
nclass<-rep(0,1000)

for (i in 1:1000){
	x<-as.numeric(duke.data[i,])

	# Recherche d'un mélange gaussien
	BIC.x<-mclustBIC(x)
	x.model<-mclustModel(x,BIC.x)
	nclass[i]<-length(x.model$parameters$mean)
}


for (i in 1:9){
	plot(density(as.numeric(duke.data[i,])), xlim=range(2,15))
	legend("topright",legend=nclass[i],bty="n")
	legend("topleft",legend=round(var(as.numeric(duke.data[i,])),3),bty="n")
	}

intercepts<-rep(0,nrow(duke.data))
for (i in 1:nrow(duke.data)){
	x<-as.numeric(duke.data[i,])
	px<-pente(x)
	intercepts[i]<-length(which(px$sign<0))
	}

i=6
x<-as.numeric(duke.data[i,])
px<-pente(x)
intercept<-which(px$sign<0)

par(mfrow=c(3,1))
hist(x, freq=F, nclass=25)
lines(density(x))
legend("topleft",legend=round(var(as.numeric(duke.data[i,])),3),bty="n")
abline(v=px$x[intercept], col="red")

plot(px$slope~px$x)
abline(h=0, col="blue")
abline(v=px$x[intercept], col="red")

plot(px$sign~px$x[-1])

px$x[intercept]


