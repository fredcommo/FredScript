vote_maj<-function(x)
{
 	n<-length(table_appr[1,]);
	res<-matrix("Indet",1,n);

	for(i in 1:n)
	{
		if(table_appr[1,i]<table_appr[2,i])
			res[,i]<-rownames(table_appr[2,]);
		{
			if(table_appr[1,i]>table_appr[2,i])
				res[,i]<-rownames(table_appr[1,]);
		}
	}
	res;
}