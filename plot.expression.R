plot.expression<-function(val,ylab="Scores",title="")

# Repr�sente les expressions par classe

#	val: Peut �tre un vecteur ou un data.frame des valeurs � repr�senter
{
	val<-as.matrix(val)

	# Compte le nombre de niveaux du crit�re 
	#lev<-length(levels(class))
	
	# Split les valeurs en fonction du crit�re
	#tmp.s<-split(val,class)

	# Partitionne la fen�tre graphique si plusieurs valeurs � repr�senter
	#p<-dim(val)[2]
	#ligne<-ceiling(p/3)
	#col<-ceiling(p/ligne)
	#par(mfrow=c(ligne,col))
	
	#if (p>1) title=colnames(val)	
	

	plot(val[2,]~val[1,],pch=19, cex=0.5,ylim=range(min(val),max(val)), main=title)
	
	p<-nrow(val)
	for (i in 3:p) points(val[i,]~val[i-1,], pch=19, cex=0.5, col=rainbow(nrow(val))[i])
}
 
