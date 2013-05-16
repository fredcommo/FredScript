# Test Affy workflow
scriptPath = "/Users/fredcommo/Documents/Projet Safir/Safir R Sources/CGHObjects/"
setwd(scriptPath)
source('SourceCodeObj.R')

# Convert a data.frame into markdown
.getWiki <- function(synId){
  synRestGET(sprintf("/entity/%s/wiki/", synId))            
}
.updateWiki <- function(wiki, synId){
  synRestPUT(sprintf("/entity/%s/wiki/%s", synId, wiki$id), wiki)
}

convert2mkd <- function(myTab, headers = TRUE){
  output <- lapply(1:nrow(myTab), function(i){paste0(paste0(as.matrix(myTab[i,]), collapse='|'),'\n')})
  output <- do.call(paste0, output)
  if (headers){
    Head <- paste0(paste0(names(Table), collapse='|'), '\n')
    output <- paste0(Head, output)
    }
  return(output)
}
ncbiLink <- function(symbol, id){
  link <- paste0('[', symbol, ']', '(http://www.ncbi.nlm.nih.gov/gene?term=', id,'[UID])')
  return(link)
}

e1  <- loadEntity('syn1864348')
source(file.path(e1$cacheDir, e1$files))
e2  <- loadEntity('syn1864353')
source(file.path(e2$cacheDir, e2$files))
e3  <- loadEntity('syn1864359')
source(file.path(e3$cacheDir, e3$files))
cgh  <- loadEntity('syn1864342')
load(file.path(cgh$cacheDir, cgh$files))
cghProfile

# create a list of symbols
geneList = c("CCND1", "ALK", "MDM2", "FRS2", "MET", "RPTOR", "ESR1", "PGR", "FGFR1", "FGFR2",
             "MYC", "FGF4", "FGF9", "EGFR", "ERBB2", "TOP2A", "IGF1", "IGF1R", "BRCA1", "BRCA2",
             "NOTCH4", "VEGFA", "PTEN", "PIK3CB", "PAK1")
# retrieve the annotations
geneTable <- geneOfInt(cghProfile, geneList)
Table <- cbind.data.frame(Symbol = paste0('**',ncbiLink(geneTable$Symbol, geneTable$entrezgene), '**'),
                          entrezgene = geneTable$entrezgene,
                          Description = geneTable$Description,
                          Chr = geneTable$Chr,
                          Cytoband = geneTable$Cytoband,
                          Log2Ratio = paste0('**',round(geneTable$Log2Ratio, 5), '**'))

#synId <- 'syn1864121'
#uri <- sprintf("/entity/%s/wiki/", synId)
#wiki <- synRestGET(uri)
synId <- 'syn1864121'
wiki <- .getWiki(synId)
mkd <- unlist(strsplit(wiki$markdown, '\n'))
mkd <- paste('${image?fileName=ExamplePlot%2Epng&align=None&scale=120}',
             '#Genes Table', convert2mkd(Table), paste0(mkd, collapse='\n'), sep ='\n')
# Replace old markdown
wiki$markdown <- mkd
wiki <- .updateWiki(wiki, synId)
#wiki <- synRestPUT(sprintf("/entity/%s/wiki/%s", synId, wiki$id), wiki)
