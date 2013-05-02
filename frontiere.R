frontiere<-function(mu1,mu2,s1,s2,p1,p2)
{
	mu1<-as.matrix(mu1,nrow=2)
	mu2<-as.matrix(mu2,nrow=2)

	solve(s2)->solve2
	solve(s1)->solve1

	l<-log(p1*sqrt(det(solve2))/p2*sqrt(det(solve1)))
	div<-t(mu2)%*%solve2%*%mu2+t(mu1)%*%solve1%*%mu1
	div<-as.numeric(div)
	cste<-(mu2-mu1)/2+(mu2-mu1)*l/div

	t(X)%*%(solve2-solve1)%*%X-(t(solve2%*%mu2)-t(solve1%*%mu1))%*%(X+cste)

	delta<-function(x)((4*x+4)^2-4*(x^2-4*x+5.562))

	x<-seq(0.126186,4.5,len=1000)
	delta(x)->dx
	b<--(4*x+4)
	y1<-(-b-sqrt(dx))/2
	y2<-(-b+sqrt(dx))/2
	lines(y1~x,col="darkgreen")
	lines(y2~x,col="darkgreen")

}

# used in plot.cloud