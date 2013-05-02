kinome.VarTest<-function(X)
{

				##### En construction ######

	# Classer X
val<-sort(X)

	# Définir l'index.0 de départ tel que index.0 = {index.valeurs| valeurs = médiane}

index.0<-which(val==median(val))

	# Remplacer les valeurs 0 par rnorm(mu=0,sd=1)

val[index.0]<-rnorm(length(index.0))
		
	# Initialise p.val
p.val<-1

	# tant que le rapport de variances n'est pas différent de 1
		# Tant que les bornes ne sont pas atteintes
	
		while(p.val>=0.05 | length(index.0)<length(val))	
		{

		# Ajouter la première valeur à gauche
		index.1<-c(min(index.0)-1,index.0)

		# tester le rapport de variances

		test<-var.test(val[index.1],val[index.0])
		p.val<-test$p.value

			# si rapport ~ 1 ajouter valeur à droite -> nouvel index de valeurs
			
				if (p.val>0.05)
				{
					index.0<-index.1
					index.1<-c(index.0,max(index.0)+1)
					test<-var.test(val[index.1],val[index.0])
					p.val<-test$p.value

					if (p.val>0.05) index.0<-index.1
				}
		if (p.val<=0.05)break()
		}
		list(min=val[min(index.0)],max=val[max(index.0)],p.value=p.val,)
}