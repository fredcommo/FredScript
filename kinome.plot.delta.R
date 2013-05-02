kinome.plot.delta<-function(X,seuil.p=0.05,title="delta")

{
n=53
p=18

# matrice des références de puits
slide.loc<-matrix("x",n,p)

	# Plage1
		plage=c(2:13)
		pas=0
		
		for (i in seq(1,15,by=2))
		{
			slide.loc[plage,18-i]<-paste(LETTERS[ceiling(i/4)+pas],seq(1,12),sep="")
			slide.loc[plage,18-(i+1)]<-paste(LETTERS[ceiling(i/4)+pas],seq(1,12),sep="")
			i=i+1
		}

		for (i in seq(4,16,by=4))
		{
			slide.loc[plage,18-(i-1)]<-paste(LETTERS[ceiling(i/4)+pas],seq(13,24),sep="")
			slide.loc[plage,18-i]<-paste(LETTERS[ceiling(i/4)+pas],seq(13,24),sep="")
			i=i+1
		}

	# Plage2
		plage=c(15:26)
		pas=4

		for (i in seq(1,15,by=2))
		{
			slide.loc[plage,18-i]<-paste(LETTERS[ceiling(i/4)+pas],seq(1,12),sep="")
			slide.loc[plage,18-(i+1)]<-paste(LETTERS[ceiling(i/4)+pas],seq(1,12),sep="")
			i=i+1
		}

		for (i in seq(4,16,by=4))
		{
			slide.loc[plage,18-(i-1)]<-paste(LETTERS[ceiling(i/4)+pas],seq(13,24),sep="")
			slide.loc[plage,18-i]<-paste(LETTERS[ceiling(i/4)+pas],seq(13,24),sep="")
			i=i+1
		}

	# Plage3
		plage=c(28:39)
		pas=8

		for (i in seq(1,15,by=2))
		{
			slide.loc[plage,18-i]<-paste(LETTERS[ceiling(i/4)+pas],seq(1,12),sep="")
			slide.loc[plage,18-(i+1)]<-paste(LETTERS[ceiling(i/4)+pas],seq(1,12),sep="")
			i=i+1
		}

		for (i in seq(4,16,by=4))
		{
			slide.loc[plage,18-(i-1)]<-paste(LETTERS[ceiling(i/4)+pas],seq(13,24),sep="")
			slide.loc[plage,18-i]<-paste(LETTERS[ceiling(i/4)+pas],seq(13,24),sep="")
			i=i+1
		}

	# Plage4
		plage=c(41:52)
		pas=12

		for (i in seq(1,15,by=2))
		{
			slide.loc[plage,18-i]<-paste(LETTERS[ceiling(i/4)+pas],seq(1,12),sep="")
			slide.loc[plage,18-(i+1)]<-paste(LETTERS[ceiling(i/4)+pas],seq(1,12),sep="")
			i=i+1
		}
		
		for (i in seq(4,16,by=4))
		{
			slide.loc[plage,18-(i-1)]<-paste(LETTERS[ceiling(i/4)+pas],seq(13,24),sep="")
			slide.loc[plage,18-i]<-paste(LETTERS[ceiling(i/4)+pas],seq(13,24),sep="")
			i=i+1
		}

# Calcul des Scores et représentation

	slide.value<-matrix(0,n,p)
	
	plot(slide.value,xlim=range(1,p),ylim=range(1:n),type="n",
			xlab="columns",ylab="rows",main=title)	

	for (i in 1:n)
		for (j in 1:p)
		{
			well<-slide.loc[i,j]
			
			if(well=="x") points(x=j,y=(54-i),pch=19,cex=1,col="lightblue")
				else 
				{
					pval<-X$p.value[which(X$well==well)]
			
					if (pval>seuil.p)
						{
						points(x=j,y=(54-i),pch=19,cex=1,col="white")
						}
					else 
					{
						score<- X$foldchange[which(X$well==well)]
						{
							if (score<0) point.color="green"
							else     point.color="red"
						}	
					points(x=j,y=(54-i),pch=19,cex=log10(abs(score))/3,col=point.color)
					}
				}
		}	
}
