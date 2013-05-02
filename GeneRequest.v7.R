
GeneRequest.v7 <- function(geneList, bySymbol = TRUE, DB = 'gene'){
  arrayInfoPath = '/Users/fredcommo/Documents/Projet Safir/Arrays Infos/'
  geneDB <- readRDS(paste0(arrayInfoPath, 'myGeneDB_2013_Mar_26.rds'))
  if(any(duplicated(geneDB$Symbol))) geneDB <- geneDB[-which(duplicated(geneDB$Symbol)),]
  if(bySymbol){
    geneDB <- geneDB[order(geneDB$Symbol),]
    output <- geneDB[is.element(geneDB$Symbol, geneList),]
    }
  else{
    geneDB <- geneDB[order(geneDB$entrezgene),]
    output <- geneDB[is.element(geneDB$entrezgene, geneList),]
    }
  return(output)
  }
