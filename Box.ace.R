Box.ace<-function(results,index)

# results : R�sultat de ACE2 (liste compl�te) 
# index : index du g�ne � repr�senter

	{
	i=index
	grp.cgh<-results$CGH.grps
	expr.values<-results$Sub.Expr[,-c(1:4)]

	grps<-which(grp.cgh[i,]!="NA")
	
	tmp<-cbind.data.frame(grp=as.character(grp.cgh[i,grps]),
			val=as.numeric(expr.values[i,grps]))

	title<-paste(results$Scores[i,1],"(",results$Scores[i,2],")")

	boxplot(tmp$val~tmp$grp,col=c("lightblue","mistyrose"),
		main=title,ylab="Log10(Expr.ratios)",xlab="CGH status")

	}
