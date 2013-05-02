HillFit<-function(x, y, graph.title=NULL){

# Calcul de l’alignement
nlm.fit<-function(bottom, top, Cycle50, slope, x){
	y.fit<-bottom + (top-bottom)^slope/(1+(Cycle50/x)^slope)
	return(y.fit)
}

# Fonction sce (somme carré résidus)
sce <- function(param, x, yobs) {
	b<-param[1]
	t<- param[2]
	C50<-param[3]
	s<-param[4]
	
	ytheo <- nlm.fit(bottom=bottom, top=t, Cycle50=C50, slope=s, x)
	return(sum((yobs - ytheo)^2))
}

Fact=0.5
X<-c(x[1]*0.5, x, x[length(x)]*4)
Y<-c(floor(y[1]), y, ceiling(y[length(y)]))

if (Y[1]<Y[2]/3) Y[1]<-Y[2]/3
bottom<-min(Y)
top<-max(Y)
Cycle50<-(min(X)+max(X))/2
slope=1

Y.init<-nlm.fit(bottom, top, Cycle50, slope, X)
plot(Y[-length(Y)]~X[-length(X)], main=graph.title)
lines(Y.init~X)

best<-nlm(f = sce, p = c(bottom, top, Cycle50, slope), x = X, yobs = Y, steptol=1e-12)
best

bottom<-best$estimate[1]
top= best$estimate[2]
Cycle50<-best$estimate[3]
slope<-best$estimate[4]

newX<-seq(min(X),max(X),length=100)
Y.best<- nlm.fit(bottom, top, Cycle50, slope, x=newX)
lines(Y.best~newX,col="blue")

Resp50<-nlm.fit(bottom,top,Cycle50,slope,Cycle50)
points(Resp50~Cycle50,col="red",pch=8,cex=2)
Resp50
#return (Resp50)
}
