kmeans.nstart<-function(x,nbcenters,iter.max=100,nstart=10)
{
	critere<-Inf;


	for(i in 1:nstart)
	{
		cur_res<-kmeans(x,nbcenters,iter.max);
		cur_critere<-sum(cur_res$withinss);
		
		if(cur_critere<critere)
		{
			critere<-cur_critere;
			res<-cur_res;
		}
	}

	res;
}
