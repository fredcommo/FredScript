Plot.means<-function(class,val)

#### En cours de construction !


mu1<- mean(bps[bps$response12==1,-c(1:2)],na.rm=T)
sd1<- sd(bps[bps$response12==1,-c(1:2)],na.rm=T)
mu2<- mean(bps[bps$response12==2,-c(1:2)],na.rm=T)
sd2<- sd(bps[bps$response12==2,-c(1:2)],na.rm=T)
n1<-dim(bps[bps$response12==1,])[1]
n2<-dim(bps[bps$response12==2,])[1]
t1<-qt(0.025,n1-1)
t2<-qt(0.02,n2-1)
IC1<-t1*sd1*sqrt(1/n1)
IC2<-t2*sd2*sqrt(1/n2)

x<-seq(1,length(mu1))
plot(mu1~I(x+0.1),xlim=range(0,length(x)+1),
ylim=range(min(mu1,mu2)-5, max(mu1,mu2)+5),pch=8,cex=1.5)
segments(x0=x+0.1,y0=mu1-IC1,x1=x+0.1,y1=mu1+IC1)
segments(x0=x-0.1,y0=mu1-IC1,x1=x+0.3,y1=mu1-IC1)
segments(x0=x-0.1,y0=mu1+IC1,x1=x+0.3,y1=mu1+IC1)

points(mu2~I(x-0.1),pch=19, cex=1.5)
segments(x0=x-0.1,y0=mu2-IC2,x1=x-0.1,y1=mu2+IC2)
segments(x0=x-0.3,y0=mu2-IC2,x1=x+0.1,y1=mu2-IC2)
segments(x0=x-0.3,y0=mu2+IC2,x1=x+0.1,y1=mu2+IC2)

legend(„topleft“,legend=c(„Non.resp“,“Resp“),pch=c(8,19),bty=“n“)

