ExtractGene <-function(GeneSet.list, base=BSets){
	
	index<- GeneSet.list$index
	output<-data.frame()
	
	for (i in 1:length(index)){
		tmp.genes<-	as.vector(base[[index[i]]]@geneIds)
		tmp.setname<- as.vector(base[[index[i]]]@setName)
		tmp<- cbind.data.frame(setname=tmp.setname, gene= tmp.genes)
		output<-rbind.data.frame(output,tmp)
		}
	output
}
