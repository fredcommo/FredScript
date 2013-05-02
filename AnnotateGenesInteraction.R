

gsearch <- function (geneSymb, database){

	# geneSymb : use official symbols. Multiple requests are accepted, e.g. "EGFR, Homo sapiens"
	# database : let's have a look at 'http://eutils.ncbi.nlm.nih.gov/entrez/query/static/eutils_help.html' for details
	# 		on available databases and other e-tools as well.
	# ! This function can return more than one Id !

	gsrch.stem <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
	gsrch.mode <- paste("db=", database, "&retmode=xml","&term=", sep = "")
	doc <- xmlTreeParse(paste(gsrch.stem, gsrch.mode, geneSymb, sep = ""), isURL = TRUE, useInternalNodes = TRUE)
	sapply(c("//Id"), xpathApply, doc = doc, fun = xmlValue)
	}



############################################
#Interaction Protéine-Protéine
############################################
get.ppiNCBI <- function(IdList) {
	
	ppi <- data.frame()
	for(i in 1:length(IdList)){
		nInt = 0
		o <- try(htmlParse(paste("http://www.ncbi.nlm.nih.gov/gene/", IdList[i], sep='')), silent = T)
		k = 1
		while(class(o)[1]=="try-error" & k<=10){
			o <- try(htmlParse(paste("http://www.ncbi.nlm.nih.gov/gene/", IdList[i], sep='')), silent = T)
			k = k + 1
			Sys.sleep(1)
			}

		# check if interaction table exists
		exist <- class(o)[1]!="try-error"
		if(exist){
			node <- node <- length(try(getNodeSet(o, "//table//th[@id='inter-prod']"), silent = T))>0
			if(node){
				p <- getNodeSet(o, "//table[@summary='Gene Interaction']")
				int <- readHTMLTable(p[[1]])
				tmp <- data.frame(egID = IdList[i], intSymbol = int$`Other Gene`)
				nInt <- nrow(tmp)
				ppi <- rbind(ppi, tmp)
				}
			# play nice! and avoid being kicked out from NCBI servers
			Sys.sleep(0.1)
			}
 
			if(dim(ppi)[1]>0){
				ppi <- unique(ppi)
				# cat(paste(dim(ppi)[1], "interactions found"), "\n")
				}
			if(nInt>0) cat(i, ":", nInt, "interactions found", "\n")
			else{
				cat(i, ":", "No interaction found", "\n")
				}
			}
		return(ppi)
}# end

require(tcltk)
require(XML)
require(org.Hs.eg.db)

glist <- annotTable$Comment.GeneName.
glist <- as.character(glist)
id <- mget(glist, org.Hs.egALIAS2EG, ifnotfound = NA)		# Gene Id
id <- unlist(id)
id <- as.numeric(as.vector(id))
nnId <- id[!is.na(id)]
# ppi <- get.ppiNCBI(id[c(2,3,5,6)])

PPI <- c()
Total = length(nnId)
pb <- tkProgressBar(title = "The R'atWork BaBar", min = 0, max = Total, width = 500)
for(i in 1:length(nnId)){
		setTkProgressBar(pb, i, label = paste("Take it easy... I'm workin' for U... :p ", round(i/Total*100, 2), "% done"))
		tmp <- get.ppiNCBI(nnId[i])
		PPI <- rbind(PPI, tmp)
		}
		close(pb)

ppi <- get.ppiNCBI(id[!is.na(id)])


PPI <- c()
for(i in 1:length(glist[1:5])){
	gname <- as.character(glist[i])
	ids <- unlist(gsearch(gname, "gene"))
	ppi <- NULL
	if(!is.null(ids) & length(ids)>0) ppi <- get.ppiNCBI(ids[1])
	PPI <- rbind(PPI, ppi)
	}

## Annotate the gene list with  homo sapiens metadata
annotate.ppi<- function(PPI){
	require(org.Hs.eg.db)
	if(!is.na(PPI$egID)){
		PPI$egSymbol <- mget(as.character(PPI$egID), envir = org.Hs.egSYMBOL, ifnotfound = NA)
		PPI$intID <- mget(as.character(PPI$intSymbol), envir = org.Hs.egSYMBOL2EG, ifnotfound = NA)
		PPI <- PPI[,c(3,2,1,4)]
		PPI
		}
	}

annotPPI <- PPI
annotPPI <- annotate.ppi(annotPPI)
annotPPI2 <- c()
for(i in 1:ncol(annotPPI)) annotPPI2 <- cbind(annotPPI2, unlist(as.character(annotPPI[,i])))
annotPPI2 <- as.data.frame(annotPPI2)
colnames(annotPPI2) <- colnames(annotPPI)

# write.table(annotPPI2, "SignaturePCA_InterAnnot.txt", sep = "\t", row.names = F)

 



