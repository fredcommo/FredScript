

library(XML)
library(annotate)


term <- "clusterin"
pmid <- esearch(term); pmid


gsearch <- function (geneId){
	gsrch.stem <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
	gsrch.mode <- "db=gene&term="
	doc <- xmlTreeParse(paste(gsrch.stem, gsrch.mode, geneId, sep=""),isURL = TRUE, useInternalNodes = TRUE)
	sapply(c("//Id"), xpathApply, doc = doc, fun = xmlValue)
	}

gsummary <- function (id){
	sum.stem <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
	sum.mode <- "db=gene&id="
	doc <- xmlTreeParse(paste(sum.stem, sum.mode, id, sep=""), isURL = TRUE, useInternalNodes = TRUE)
	sapply(c("//Item"), xpathApply, doc = doc, fun=xmlValue)
	}

esearch <- function (term){
	srch.stem <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
	srch.mode <- "db=pubmed&retmax=1000&retmode=xml&term="
	doc <-xmlTreeParse(paste(srch.stem,srch.mode,term,sep=""), isURL = TRUE, useInternalNodes = TRUE)
	sapply(c("//Id"), xpathApply, doc = doc, fun = xmlValue)
	}


geneId <- 'Hs.635482'	# clusterin CLU
gene <- gsearch(geneId); gene
id <- unlist(as.numeric(gene))

id = 2257			# FGF12 entrezGeneId = 2257
test <- gsummary(id); test
gene.symb <- as.character(as.table(test)[1])
gene.name <- as.character(as.table(test)[2])

pmid <- esearch(paste("(", gene.symb, "OR", gene.name, ")", "AND docetaxel")); pmid
pmid <- unlist(pmid)

x<- pubmed(pmid)
a <- xmlRoot(x)
numAbst <- length(xmlChildren(a)); numAbst

arts <- vector("list", length = numAbst)
abstrs <- rep(NA, numAbst)

	for(i in 1:numAbst){
		arts[[i]] <- buildPubMedAbst(a[[i]])
		abstrs[i] <- abstText(arts[[i]])
	}


absts <- list()
for (i in 1:numAbst) {
	absts[[i]] <- buildPubMedAbst(a[[i]])
	}

## file for the output

output.title <- paste(gene.name, "PubMed.html")
pmAbst2HTML(absts,filename=output.title)


# bon format
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmax=1000&retmode=xml&term=clusterin
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=HS.436657
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nucleotide&term=NM_203339
http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=	# output gsearch() = 1191