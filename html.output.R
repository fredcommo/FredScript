html.output<-function(data, Chip="hgu133a", n.rows=nrow(data), html=TRUE, title="")

# For html output

# data	: data frame to annotate. Columns have to be named as data$t.value, data$adj.p.value
# Chip	: The Chip database to consider
# n.rows	: The number of results to consider
# title	: The title for the table

{
	library(annaffy)
		
		data <- data[1:n.rows,]
		probeId<-rownames(data)
		probeId<-as.character(probeId)
	
		value.table<-aafTable("t.value" = signif(data$t.value, 3),
						"Adj.p.value" = signif(data$adj.p.value, 3),
						signed=TRUE)

		# 	1: Probe	2: Symbol	3: Description	4: Chromosome	5: Chromosome Location
		#	6: GenBank	7: Gene	8: Cytoband		9: UniGene		10: PubMed             
		#	11: Gene Ontology		12: Pathway
 
		ann.col<-aaf.handler()[c(1:3,7,8,10:12)]

		ann.table<-aafTableAnn(probeId, Chip, ann.col)
		ann.table<-merge(ann.table, value.table)
	
		saveHTML(ann.table, paste(title, ".html", sep=""), title=title)

}

