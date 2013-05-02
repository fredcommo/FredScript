ROC.density<-function(class, val, cut, add=F, title=NULL)

{
class <- as.factor(class)
lev<-levels(class)
x1<-density(val[class==lev[1]],na.rm=T)$x
y1<-density(val[class==lev[1]],na.rm=T)$y

x2<-density(val[class==lev[2]],na.rm=T)$x
y2<-density(val[class==lev[2]],na.rm=T)$y

range.x<-c(-max(y1,y2),max(y1,y2))
range.y<-c(min(x1,x2),max(x1,x2))

y1<-(-y1)

if(!add) plot(x1~y1,type="l",xlim=range.x,ylim=range.y,col="blue",xlab="Density",ylab="values",main=title)
lines(x2~y2,col= "red")
abline(v=0,lty=3)
abline(h=cut,lty=1,col="blue")

legend(x=range.x[1],y=cut,legend=cut,text.col="blue",cex=.75,bty="n",y.intersp=-2.5)
legend("bottomleft",legend=lev[1],text.col= "blue",bty= "n")
legend("bottomright",legend=lev[2],text.col= "red",bty= "n")
}
