require(annotate)
require(XML)
require(biomaRt)
require(foreach)
require(iterators)

# url example : http://www.genome.jp/kegg-bin/get_htext?htext=br08303_target.keg&query=EGFR

# build a KeggToDrug repository
repofun <- function(ids, ...)
  paste0("http://www.genome.jp/kegg-bin/get_htext?htext=br08303_target.keg&query=", ids)
setRepository("keggD", repofun)

itemList = c('NomenclatureSymbol', 
             'NomenclatureStatus', 'OtherAliases', 'Chromosome',
             'MapLocation', 'Orgname')

# gSearch: return gene Ids
gSearch <- function (geneSymb, database){  
  '
	Called by geneRequest.v7()
  '
  # ! This function can return more than one Id ! 
  gsrch.stem <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
  gsrch.mode <- paste0("db=", database, "&retmode=xml","&term=")
  URL <- paste0(gsrch.stem, gsrch.mode, geneSymb)
  doc <- xmlTreeParse(URL, isURL = TRUE, useInternalNodes = TRUE)
  sapply(c("//Id"), xpathApply, doc = doc, fun = xmlValue)
}

getItem <- function(doc, item){
  '
	Called by gSummary()
  '
	expr = paste0('<Item Name=\"', item, '\"')
  if(any(grepl(expr, doc))){
    String = doc[grep(expr, doc)]
    r = regexec('>(.*?)<', String)
    return(unlist(regmatches(String, r))[2])
  }
  return(NA)
}

gSummary <- function(id, database){
  '
	Called by geneRequest.v7()
  '
	sum.stem <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
  sum.mode <- paste0("db=", database, "&id=")
  urlFile = url(paste0(sum.stem, sum.mode, id), 'r')
  doc = readLines(urlFile)
  close(urlFile)
  geneInfo = c()
  foreach(item = iter(itemList)) %do% {
    geneInfo = cbind(geneInfo, getItem(doc, item))
    colnames(geneInfo)[length(geneInfo)] = item
  }
  return(as.data.frame(geneInfo))
}

geneRequest.v7 <- function(geneList, database, verbose = TRUE){
  output = data.frame()
  foreach(gene = iter(geneList)) %do% {
    notFound = cbind.data.frame(gene, t(rep(NA, length(itemList))))
    colnames(notFound) = c(itemList, 'entrezgene')
    id = gSearch(paste0(gene, '[symbol]+homo+sapiens[Organism]'), database)
    id = unlist(id)
    if(length(id) == 0){
      if(verbose) cat('\n', gene, '\t*** not found ***')
      output = rbind.data.frame(output, notFound)
    }
    else{
      # should have NomenclatureStatus = 'Official'
      Official = FALSE
      k = 1
      while (!Official & k <= length(id)){
        tmp = gSummary(paste0(id[k], '%5BUID%5D'), database)
        tmp = cbind.data.frame(tmp, entrezgene = id[k])
        Official <- tmp$NomenclatureStatus == 'Official'
        k = k + 1
      }
      if(!Official) tmp = notFound
      output = rbind.data.frame(output, tmp)
      if(verbose) cat('\n', gene, 'found')
    }
  }
  rownames(output) = seq(1, nrow(output))
  return(output)
}

KeggToDrug <- function(ID){
  srch.stem <- "http://www.genome.jp/kegg-bin/get_htext?"
  srch.mode1 <- "htext=br08303_target.keg&query="
  srch.mode2 <- "&htext=br08303_target.keg"
  loc <- paste(srch.stem, srch.mode1, ID, srch.mode2, sep="")
  doc <- htmlParse(loc)
  tabNodes <- getNodeSet(doc, "//table")
  tabHtml <- readHTMLTable(tabNodes[[3]])
  drug.list <- as.character(tabHtml[-1,2])
  ndrug <- length(which(!is.na(drug.list) & drug.list != ""))
  return(cbind.data.frame(ndrug = ndrug, htmlLink = loc))
}

# Load full gene symbols HUGO list


human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
output = getBM(attributes=c('hgnc_symbol', 'description', 'chromosome_name', 'band', 'start_position','end_position', 'entrezgene'),
              mart = human)
output$description =  gsub(' \\[.*\\]', '', output$description)
output = as.data.frame(output)
#nrow(output)

keggList <- data.frame()
foreach(i = 1:nrow(output)) %do%{
  ID = as.character(output$entrezgene[i])
  Symbol = output$hgnc_symbol[i]
  if(is.na(ID) & Symbol != ''){
    webReq = geneRequest.v7(Symbol, 'gene', verbose = FALSE)
    ID = as.character(webReq$entrezgene)
    }
  if(!is.na(ID) & nchar(ID)>3)
    keggList = rbind.data.frame(keggList, cbind.data.frame(geneId = ID, Symbol = Symbol, KeggToDrug(ID)))
  else if(Symbol!='')
    keggList = rbind.data.frame(keggList, cbind.data.frame(geneId = ID, Symbol = Symbol, KeggToDrug(Symbol)))
  if(i%%100 == 0) cat(i, 'of', nrow(output))
  if(i%%1000 == 0) cat('\n')
}



write.table(keggList, "Kegg_Drugs_Table.txt", sep = "\t", row.names = F)

