ACE2<-function(ACE.res,expr,cgh)

# ACE.results : fichier issu de l'analyse ACE-it (edit results)
# expr : Fichier GE utilisé pour ACE-it (donc de même forme)
# cgh : Fichier des données CGH utilisé pour ACE-it (donc de même forme)

# Contenu des fichiers GE et CGH :
#		    ID	chromosome	   start	   end     Values	... ->
#	A_14_P112718	         1	  554268	554327	    0
#	A_14_P201353	         1	  749622	749681	    0


{
	best<-ACE.res
	best$RefCGH->refcgh		# ref des sondes CGH
	best$RefGE->refexpr		# ref des sondes GE

	index.expr<-which(expr$Id %in% refexpr)
	sub.expr<-expr[index.expr,]	# sous-tableau des sondes GE-ACEit

	subcgh<-matrix(0,dim(best)[1],dim(cgh[,-c(1:4)])[2])
	
	for (i in 1:length(refcgh)){
		probe<-as.character(refcgh[i])
		subcgh[i,]<-as.numeric(cgh[which(cgh[,1]==probe),-c(1:4)])
		}

	subcgh<-as.data.frame(subcgh)
	colnames(subcgh)<-colnames(cgh[,-c(1:4)])
	grp.cgh<-subcgh

	for (i in 1:dim(subcgh)[1])
		{
		if (best$Comparison[i]=="l_vs_n")
			{
			grp.cgh[i,]<-ifelse(subcgh[i,]<0,"Loss",NA)
			grp.cgh[i,]<-ifelse(subcgh[i,]==0,"N",grp.cgh[i,])
			}
		else
			{ 
			grp.cgh[i,]<-ifelse(subcgh[i,]>0,"Gain",NA)
			grp.cgh[i,]<-ifelse(subcgh[i,]==0,"N",grp.cgh[i,])
			}
		}

	library(multtest)
	expr.values<-sub.expr[,-c(1:4)]
	res<-data.frame(index=0,teststat=0,rawp=0,adjp=0)

	for (i in 1:dim(grp.cgh)[1])
		{
		index<-which(grp.cgh[i,]!="NA")
		grp<-as.character(grp.cgh[i,index])
		grp<-c(0,1)[factor(grp)]
		res[i,]<-mt.maxT(expr.values[i,index],classlabel=grp,nonpara="y")
		}

	res[,4]<-p.adjust(res[,3],method="BH")
	res<-cbind(Gene=best[,3],Band=best[,2],res)
	res$index<-seq(1,dim(res)[1])

	list(Scores=res,Sub.Expr=sub.expr,Sub.CGH=subcgh,CGH.grps=grp.cgh)
}

