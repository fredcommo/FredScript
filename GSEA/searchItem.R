searchItem<-function(Search="",base=BNames){
	output<-data.frame()
	p<-dim(base)[1]
	#Search<- toupper(Search)
	for (i in 1:p)
{
	items<-unlist(strsplit(base$setname[i],"_"))
	if (Search %in% items) output<-rbind.data.frame(output,base[i,])
	}
	output
}

