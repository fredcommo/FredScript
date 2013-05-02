conf<-function(...)
{
	tab<-table(...);
	nb<-sum(tab);
	
	erreur<-round(100*(sum(tab)-sum(diag(tab)))/sum(tab),2);

	cat("Taux d'erreur= ",erreur,"%\n");
	erreur;
}