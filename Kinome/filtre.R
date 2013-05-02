filtre<-function(Results,min.p,min.fc,min.int)
{
	value1<-Results$estim1
	value2<-Results$estim2
	fc<-Results$fc
	p.value<-Results$p.value

	index.up<-which(fc>=min.fc & p.value<=min.p & (value1>min.int | value2>min.int))
	index.low<-which(fc<=(-min.fc) & p.value<=min.p & (value1>min.int | value2>min.int))

	list(Results, up=Results[index.up,], low=Results[index.low,])
}

