CVqda<-function(train,V,...)
{
	library(MASS);
	rand<-sample(V,dim(train)[1],replace=T);
	res<-train$letter;
	for(i in sort(unique(rand)))
	{
		appr<-qda(letter~.,train[rand!=i,]);
		res[rand==i]<-predict(appr,train[rand==i,])$class;
	}
	res;
}