kinome.rename<-function(X,File.infos=infos)

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

	# Ajoute les données infos correspondantes

		infos<-File.infos
		n<-dim(X)[1]
		res<-cbind.data.frame(X[1,1:2],infos[1,3:11])
		res[1,]<-NA
		for (i in 1:n)
		{
			well<-X$well[i]
			{	if (well!="Blank")
				{
					well.value<-X[i,]
					well.info<-infos[which(infos$well==well),3:11]
					res[i,]<-cbind.data.frame(well.value,well.info)
				}
				else res[i,1:2]<-X[i,]
			}
		}
		X<-res
		X
}


