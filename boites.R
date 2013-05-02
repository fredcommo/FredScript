boites<-function(x,...)
{
	data<-x;
	n<-dim(data)[2];

	boxplot(data,border="white");

	for (i in 1:n)
	{
		m<-median(data[,i]);
		
			if (m>0) 
				{
				boxplot(data[,i],names=colnames(data)[i],border="green",add=T,at=i)
				text(i,(m+1/sample(1:10,1)),labels=dimnames(band1032)[[2]][i]);
				}

			if (m<0) 
				{
				boxplot(data[,i],names=colnames(data)[i],border="red",add=T,at=i)
				text(i,(m-1/sample(1:10,1)),labels=dimnames(band1032)[[2]][i]);
				}

		if (m==0) boxplot(data[,i],names=colnames(data)[i],add=T,at=i);
		
		

	}
}
