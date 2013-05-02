microarray.formate<-function(X,Values=NA)		#F635Mean_B635

# Formatage des données, récupération des réplicats. Pour fichiers issus des scanners.


{
	if (is.na(Values[1])) stop ("Aucune valeur sélectionnée pour la mise en forme !")

	col.index<-which(colnames(lame1) %in% Values)
	p<-length(Values)
	
	well.Id<-paste(LETTERS[1],seq(1,24),sep="")

	for (i in 2:16)
		well.Id<-c(well.Id,paste(LETTERS[i],seq(1,24),sep=""))

	output<-data.frame(well=well.Id,kinase=0,substrate=0)
	
	for (i in 1:p)
		output<-cbind.data.frame(output,name1=0,name2=0)	

	colnames(output)[4:(3+p*2)]<-paste(rep(Values,each=2),rep(1:2),sep=".")

	nline<-dim(output)[1]
	
	for (i in 1:nline)
	{
		tmp<-X[X$Spot.Coord==well.Id[i],]
		output$kinase[i]<-as.character(tmp$ID)[1]
		output$substrate[i]<-as.character(tmp$Substrate)[1]
		output[i,4:(3+p*2)]<-as.numeric(unlist(tmp[,col.index]))	
	}
	output
}
