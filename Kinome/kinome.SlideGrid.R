kinome.SlideGrid<-function(n=53,p=18)

# matrice des références de puits

{

slide.loc<-matrix("x",n,p)

	# Plage1
		plage=c(2:13)
		pas=0
		
		for (i in seq(1,15,by=2))
		{
			slide.loc[plage,p-i]<-paste(LETTERS[ceiling(i/4)+pas],seq(1,12),sep="")
			slide.loc[plage,p-(i+1)]<-paste(LETTERS[ceiling(i/4)+pas],seq(1,12),sep="")
			i=i+1
		}

		for (i in seq(4,16,by=4))
		{
			slide.loc[plage,p-(i-1)]<-paste(LETTERS[ceiling(i/4)+pas],seq(13,24),sep="")
			slide.loc[plage,p-i]<-paste(LETTERS[ceiling(i/4)+pas],seq(13,24),sep="")
			i=i+1
		}

	# Plage2
		plage=c(15:26)
		pas=4

		for (i in seq(1,15,by=2))
		{
			slide.loc[plage,p-i]<-paste(LETTERS[ceiling(i/4)+pas],seq(1,12),sep="")
			slide.loc[plage,p-(i+1)]<-paste(LETTERS[ceiling(i/4)+pas],seq(1,12),sep="")
			i=i+1
		}

		for (i in seq(4,16,by=4))
		{
			slide.loc[plage,p-(i-1)]<-paste(LETTERS[ceiling(i/4)+pas],seq(13,24),sep="")
			slide.loc[plage,p-i]<-paste(LETTERS[ceiling(i/4)+pas],seq(13,24),sep="")
			i=i+1
		}

	# Plage3
		plage=c(28:39)
		pas=8

		for (i in seq(1,15,by=2))
		{
			slide.loc[plage,p-i]<-paste(LETTERS[ceiling(i/4)+pas],seq(1,12),sep="")
			slide.loc[plage,p-(i+1)]<-paste(LETTERS[ceiling(i/4)+pas],seq(1,12),sep="")
			i=i+1
		}

		for (i in seq(4,16,by=4))
		{
			slide.loc[plage,p-(i-1)]<-paste(LETTERS[ceiling(i/4)+pas],seq(13,24),sep="")
			slide.loc[plage,p-i]<-paste(LETTERS[ceiling(i/4)+pas],seq(13,24),sep="")
			i=i+1
		}

	# Plage4
		plage=c(41:52)
		pas=12

		for (i in seq(1,15,by=2))
		{
			slide.loc[plage,p-i]<-paste(LETTERS[ceiling(i/4)+pas],seq(1,12),sep="")
			slide.loc[plage,p-(i+1)]<-paste(LETTERS[ceiling(i/4)+pas],seq(1,12),sep="")
			i=i+1
		}
		
		for (i in seq(4,16,by=4))
		{
			slide.loc[plage,p-(i-1)]<-paste(LETTERS[ceiling(i/4)+pas],seq(13,24),sep="")
			slide.loc[plage,p-i]<-paste(LETTERS[ceiling(i/4)+pas],seq(13,24),sep="")
			i=i+1
		}
	slide.loc
}
