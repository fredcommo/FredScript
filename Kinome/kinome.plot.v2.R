kinome.plot.v2<-function(X,title="", flags=TRUE, suppr=TRUE, coef=1000,log.base=10)

# X			: file containing values (file.QC$Values or file.norm) and QC status.
# title		: title for the graph
# flags		: should divergent duplicates have to be suppressed. Default is TRUE.
# suppr		: should duplicates for which at least one value is lower than bg have to be suppressed. Default is TRUE.
# coef, log.base 	: graphical parameters for the sizes of spots.

{
	input<-X
	n<-dim(input)[1]
	
	# Suppression des valeurs négatives 
	input$value1<-ifelse(input$value1<0, NA, input$value1)
	input$value2<-ifelse(input$value2<0, NA, input$value2)

	if(suppr)
	{
		input$value1<-ifelse(input$QC=="suppr",NA,input$value1)
		input$value2<-ifelse(input$QC=="suppr",NA,input$value2)
	}
	

	# Si filtrage des flags => valeur du spot = NA (non calculé)
	if (flags)
	{		
		input$value1<-ifelse(input$QC=="flag", NA, input$value1)
		input$value2<-ifelse(input$QC=="flag", NA, input$value2)
	}

	# matrice des positions
	source("D:\\Stats\\Doc R\\Scripts R\\Kinome\\kinome.SlideGrid.R")
	slide.pos<-kinome.SlideGrid()
	
	# matrice des valeurs
	slide.val<-matrix(NA,53,18)	#coef <-> NA

	for (i in seq(2,16,by=2))
		for (j in 1:53)
		{
			well<-slide.pos[j,i]
			if(well!="x")
				{
				slide.val[j,i]<-input$value1[which(input$well==well)]
				slide.val[j,i+1]<-input$value2[which(input$well==well)]
				}
		}

	n<-dim(slide.val)[1]	# rows = y axis
	p<-dim(slide.val)[2]	# columns = x axis

	val.max<-max(slide.val,na.rm=TRUE)

	plot(1, 1, xlim=range(1,p), ylim=range(1:n),type="n", xlab="columns", ylab="rows", main=title)	
	
	#version4 (voir programme R)
	for (i in 1:p)
		for (j in 1:n)
		{	
			val<-log(slide.val[j,i]/coef,base=log.base)
			if (!is.na(val) & val < 0) val=0.1
			points(x=i,y=(54-j),pch=19,cex=val)
		}
	
	# Repères
	for (j in 1:53)
	{
		points(x=1,y=j,pch=19,cex=0.5,col="lightblue")
		points(x=18,y=j,pch=19,cex=0.5,col="lightblue")
	}
	for (j in c(1,14,27,40,53))
		for (i in seq(1,18))
			points(x=i,y=j,pch=19,cex=0.5,col="lightblue")
	
	# fichier de sortie des valeurs
	slide.val
}
