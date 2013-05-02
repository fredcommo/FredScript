kinome.rename1<-function(X)

{
	source("D:\\Stats\\Doc R\\Scripts R\\Kinome\\kinome.SlideGrid.R")
	grid<-kinome.SlideGrid()
	
	# renomme les positions
		X$well<-as.vector(X$well)

		n<-dim(grid)[1]
		p<-dim(grid)[2]

		j=1
		for (i in p:1)
		{
			X$well[j:(j+n-1)]<-grid[,i]
			j=j+n
		}
		X$well[which(X$well=="x")]<-"Blank"

		X
}


