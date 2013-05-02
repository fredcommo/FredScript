

pm.srch<- function (){
srch.stem <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term="
query <-as.character(scan(file="",what="character"))
doc <-xmlTreeParse(paste(srch.stem,query,sep=""),isURL = TRUE,
useInternalNodes = TRUE)
sapply(c("//Id"), xpathApply, doc = doc, fun = xmlValue)
}

doc <- pm.srch()
"hif3a"
 
#####################################

library(XML)
library(annotate)

esearch <- function (term){
	srch.stem <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
	srch.mode <- "db=pubmed&retmax=1000&retmode=xml&term="
	doc <-xmlTreeParse(paste(srch.stem,srch.mode,term,sep=""),isURL = TRUE, useInternalNodes = TRUE)
	sapply(c("//Id"), xpathApply, doc = doc, fun = xmlValue)
	}

term <- "clu AND docetaxel"
pmid <- esearch(term); dim(pmid)
pmid <- unlist(pmid)

x<- pubmed(pmid)
a <- xmlRoot(x)
numAbs <- length(xmlChildren(a)); numAbs

arts <- vector("list", length = numAbs)
abstrs <- rep(NA, numAbs)

for(i in 1:numAbs){
	arts[[i]] <- buildPubMedAbst(a[[i]])
	abstrs[i] <- abstText(arts[[i]])
	}

# A finir pour sorties html choix multiples/complémentaires

found <- grep("docetaxel", abstrs)
goodAbstrs <- arts[found]
length(goodAbstrs)

test <- genbank()
ptest <- parseDTD(test)
xmlContainsEntity(PubMedId, test)
xmlContainsElement(factor(PubMedId), test)

Gene-ref_desc
Gene-ref_locus

