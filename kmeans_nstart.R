## Prend data et un nombre de centres en arguments
## Calcule nfois (nstart aléatoires) l'inertie intra-classes ($withinss) et retient la meilleure partition.


kmeans_nstart<-function(x,nbcenters,iter.max=10000,nstart=10)
{
        critere<-Inf;		# initialise à la valeur max


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