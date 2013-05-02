plot.express<-function(class, val, Ylab = "Scores", Test = TRUE, test.type = "student", title = "")

# Représente les expressions par classe
#	class: vecteur factor à discriminer
#	val: Peut être un vecteur ou un data.frame des valeurs à représenter

{
	class <- as.factor(class)
	val <- as.data.frame(val)

	# Compte le nombre de niveaux du critère 
	lev <- length(levels(class))
	
	# Split les valeurs en fonction du critère
	tmp.s <- split(val,class)

	# Partitionne la fenêtre graphique si plusieurs valeurs à représenter
	p <- dim(val)[2]
	if(p>1){
		ligne <- ceiling(p/3)
		col <- ceiling(p/ligne)
		par(mfrow = c(ligne,col))
		}

	if (p>1) title = colnames(val)	
	
	# de 1 à p-ème valeur	
	palette <- c("goldenrod3", "royalblue1", "purple3", "indianred3")
	for (i in 1:p)
	{
		y <- sort(tmp.s[[1]][,i])
		maxi <- max(sort(val[,i]))
		mini <- min(sort(val[,i]))
		n <- length(y)
		x <- seq(0,1,by=1/(n-1))
		plot(y~x, type="n", xlab="Fn(x)", ylim=range(mini,maxi), ylab = Ylab, main=title[i])

		for (l in 1:lev)	
		{
			y <- sort(tmp.s[[l]][,i])
			n <- length(y)
			x <- seq(0,1,by=1/(n-1))
			lines(y~x, col=palette[l], lwd = 3)
		}

		X <- cbind.data.frame(class,val)
		tab <- table(X$class)
		legend("topleft", legend=paste(levels(class)," (n=",as.vector(tab),")",sep=""),lwd = 3, col=palette[1:lev], cex = 1.25, bty="n")

		if(Test){		
			if (lev<3) 
				{
				n1 <- length(val[class==lev[1],i])
				n2 <- length(val[class==lev[2],i])

				if (n1==n2) paired.test = TRUE
				else paired.test = FALSE
				
				if(test.type=="student") {
					p.val <- t.test(val[,i]~class, paired = paired.test)$p.value
					legend("bottomright", legend=c("t-test", paste("p = ",signif(p.val,3),sep="")), cex = 1, text.col="darkblue", bty="n")
					}
				if(test.type=="wilcox"){
					p.val <- wilcox.test(val[,i]~class, paired = paired.test)$p.value
					legend("bottomright", legend=c("Wilcoxon-test", paste("p = ",signif(p.val,3),sep="")), cex = 1, text.col="darkblue", bty="n")
					}
				
				} 	
			else 
				{
				p.val <- kruskal.test(val[,i]~class)$p.value
				legend("bottomright", legend = c("Kruskal.test", paste("p=",signif(p.val,3), sep="")), cex = 1, text.col = "darkblue", bty = "n")
				}
		}
	}
}
 
