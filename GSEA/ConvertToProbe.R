ConvertToProbe<-function(gene.list, Symb, Split=F){

	setnames<- levels(gene.list$setname)
	n.list<-length(setnames)
	probe.list<-data.frame()

	for (i in 1:n.list){
		tmp.list<-gene.list[which(gene.list$setname==setnames[i]),]
		p<- nrow(tmp.list)

		for (j in 1:p){
			set.name<- setnames[i]
			gene.symb<-tmp.list$gene[j]
	probes<-as.vector(names(which(Symb==gene.symb)))
	if(length(probes)==0) probes<-NA

	tmp<- cbind.data.frame(setname=set.name, gene= gene.symb,
					probe=probes)
			probe.list<-rbind.data.frame(probe.list,tmp)
		}
	}
	if (Split) probe.list<-split(probe.list,probe.list$setname)

	probe.list
}
