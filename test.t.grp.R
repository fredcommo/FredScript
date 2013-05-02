test.t.grp<-function(X, p=0.05)

	# test.t.grp(log2,p=0.05)calcule un test t et renvoie les séquences significatives
	# log2 : un sous-groupe de log2.ratio sélectionné 
	#	selon un index which(phenotype(aCGH.obj)$critère==critère)
	# p : p-value seuil. par défaut p=0.05

	{
		grp<-X
		n<-dim(grp)[1]
		p.val<-matrix(0,n,2)
		rownames(p.val)<-rownames(grp)
		colnames(p.val)<-c("index","p.value")

		p.val[,1]<-seq(1:n)
		mu<-mean(as.matrix(grp))

		for (i in 1:n)
		{
		p.val[i,2]<-t.test(grp[i,],mu=mu,equal.var=F,conf.level=(1-p))$p.value
		}
		
		p.val<-as.data.frame(p.val)
		p.val[p.val$p.value<=p,]
	}

	# Valeurs
	# la fonction renvoie les séquences satisfaisant le seuil de p-value
	