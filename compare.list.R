compare.list<-function(list1,list2)

# Compare 2 listes et 
#	renvoie un tableau avec le nombre d'éléments communs et spécifiques de chaque liste

{
	communs<-which(list1 %in% list2)
	if (length(communs)!=0)
	{
		specif.list1<-list1[-communs]
		specif.list2<-list2[-communs]
	}
	else
	{
		specif.list1<-list1
		specif.list2<-list2
	}
	col1<-paste("specif.",substitute(list1),sep="")
	col2<-paste("specif.",substitute(list2),sep="")

	tab.res<-data.frame(length(communs), length(specif.list1), length(specif.list2))
	colnames(tab.res)<-c("communs",col1,col2)
	tab.res
	
	# list(col1=length(list1),col2=length(list2),tab.res)
}