
# Recherche les GO.ref à partir d'un terme de recherche

library(annaffy)
GOTerm2Tag<-function(term){
	GTL<-eapply(GOTERM, function(x){
		agrep(term, x@Term, value=TRUE)})

	GL<-sapply(GTL, length)
	names(GTL[GL>0])
}

# exemple :
# Repair<-GOTerm2Tag ("repair")

# [1] "GO:0051103" "GO:0045020" "GO:0045021" "GO:0042275" "GO:0042276" "GO:0045738" "GO:0045739"
# [8] "GO:0010213" "GO:0000711" "GO:0000731" "GO:0000734" "GO:0043504" "GO:0006281" "GO:0006282"
# [15] "GO:0046787"

as.character(mget(Repair,hgu133aGO2PROBE))


	library(hgu133a)
	all.probes<-rownames(all)

	# obtenir la liste des recherches possibles : ls("package:hgu133a")

	GO.terms<-unlist(mget(all.probes,hgu133a2GO))

	data.frame(get("GO:0016740",hgu133a2GO2PROBE))

 	good.terms<-GO.terms[which(GO.terms %in% Repair)]
  	gene.list<-rep(NA,length(good.terms))
	for (i in 1:length(good.terms))
		gene.list[i]<-strsplit(names(good.terms)[i],".GO:*")[[1]][1]
	gene.list<-as.data.frame(gene.list)

	sub<-data[which(rownames(data) %in% gene.list),]

	correl2<-rep(0,nrow(sub))
	for (i in 1:nrow(sub))
		correl2[i]<-cor.test(as.numeric(ercc1[,2]),as.numeric(sub[i,]),method="spearman")$estimate

	correl2<-cbind.data.frame(index=seq(1,length(correl1)),correl1)
	correl2<-correl1[order(abs(correl1[,2]),decreasing=T),]

		

	library(mclust)
	n.gene<-nrow(data)
	n.mode<-data.frame(prode=rownames(data), n=0, p.min=0)

	for (i in 1:n.gene)
	{
		x<-data[i,]
		BIC.x<-mclustBIC(x)
		x.model<-mclustModel(x,BIC.x)
		n.mode$n[i]<-length(x.model$parameters$mean)
		if (n.mode$n[i]==1) n.mode$p.min[i]<-1
		else n.mode$p.min[i]<-min(x.model$parameters$pro)
	}

	par(mfrow=c(1,2))
	mclust1Dplot(as.numeric(x), parameters = x.model$parameters, z = x.model$z,
			what = "density", identify = TRUE)
	mclust1Dplot(x, z = x.model$z, parameters = x.model$parameters,
			what = "uncertainty", identify = TRUE)


