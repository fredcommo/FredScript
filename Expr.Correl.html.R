Expr.Correl.html<-function(x, data, method="spearman", Chip="hgu133a", n.rows=50, html=TRUE, title="")

# For html output of correlation results

# x 		: reference (expression to be compared to)
# data	: data frame of all the values to compare to the reference
# method	: method used for the correlation. Default = "spearman"
# Chip	: The Chip database to consider
# n.rows	: The number of results to consider
# title	: The title for the table

{
	library(annaffy)

	output<-data.frame(Probe=rownames(data),estimate=0,p.value=0)
	x<-as.numeric(x)+rnorm(length(x),mean=0,sd=1e-10)	
		
		for (i in 1:nrow(data)){	
			y<-as.numeric(data[i,])+rnorm(length(x),mean=0,sd=1e-10)

			test<-cor.test(x,y,method="spearman",exact=TRUE)
			output$estimate[i]<-test$estimate
			output$p.value[i]<-test$p.value
		}

	output<-output[order(abs(output$estimate),decreasing=TRUE),]

	if (html){

		probeId<-output$Probe[1:n.rows]
		probeId<-as.character(probeId)
	
		values<-round(output$estimate[1:n.rows],4)
		value.table<-aafTable("Spearman.rho"=values, signed=TRUE)

		# 	1: Probe	2: Symbol	3: Description	4: Chromosome	5: Chromosome Location
		#	6: GenBank	7: Gene	8: Cytoband		9: UniGene		10: PubMed             
		#	11: Gene Ontology		12: Pathway
 
		ann.col<-aaf.handler()[c(1:3,7,8,10:12)]

		ann.table<-aafTableAnn(probeId, Chip, ann.col)
		ann.table<-merge(ann.table, value.table)
	
		saveHTML(ann.table, paste(title, ".html", sep=""), title=title)

	}

	output
}



