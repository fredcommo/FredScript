#################################################################################################

# Required packages
require(XML)      # used for xmlTreParse(xmlDoc)
require(annotate) # used for pubmed(pubmed_ids)
require(tcltk2)   # for progressBars. ######### Does not work in RStudio #########
require(R2HTML)   # to manipulate html docs

GeneRequest.v6 <- function(geneList, bySymb = TRUE, DB = "gene", verbose = TRUE){
  geneList <- as.character(geneList)
  Total = length(geneList)
  # pb <- tkProgressBar(title = "The R'atWork BaBar", min = 0, max = Total, width = 500)  
  # cum.len <- cumsum(hg19.info$length)
  
  gene.annot <- c()
  for(i in 1:Total){    
    Sys.sleep(0.1)
    # launch & increment the pBar
    # setTkProgressBar(pb, i, label = paste("Take it easy... I'm workin' for U... :p ", round(i/Total*100, 2), "% done"))    
    Symb <- "Not found"
    Name <- Org <- Chr <- Cytoband <- Chr.start <- Chr.end <- RefSeq <- Id <- NA			# in case of error only !
    ids <- geneList[i]
    
    if(bySymb){
      gname <- geneList[i]
      ids <- unlist(gSearch(paste(gname, "[symbol] homo sapiens", sep = ""), database = DB))					#homo sapiens
      if(is.null(ids)) cat("\n ***", gname, "not found: Stop kidding *** !\n\n")
      }
    if(!is.null(ids)){ 
      if(bySymb & verbose) cat(gname, "found:", length(ids), "ids.\n")
      
      j = 1
      Id <- ids[j]
      gsum <- unlist(gSummary(paste(Id, "+homo+sapiens", sep = ""), database = DB))					#homo sapiens
      foundName <- gsum[1]
      Org <- gsum[3]
      status <- gsum[2]
      verifId <- as.numeric(gsum[5])
      Chr <- gsum[6]
      
      if(length(ids)>1){
        while ((foundName != gname| Org != "Homo sapiens" | status=="reserved" | verifId != 0 | Chr =="") & j < length(ids)){
          j = j + 1
          Id <- ids[j]
          gsum <- unlist(gSummary(paste(Id, "homo sapiens"), database = DB))			#homo sapiens
          foundName <- gsum[1]
          status <- gsum[2]
          Org <- gsum[3]
          verifId <- as.numeric(gsum[5])
          Chr <- gsum[6]
          }
        # cat(gname, "related to Homo sapiens found at occurence:", j, "\n")
        }
      Symb <- gsum[1]
      if(!bySymb){
        gname <- Symb
        cat(as.character(Id), "found as", gname, "\n")
        }
      Name <- gsum[2]
      Chr <- gsum[6]
      Cytoband <- gsum[8]
      if(Chr == "") Chr <- NA
      if(!is.na(Chr) & Chr == "X") Chr <- 23
      if(!is.na(Chr) & Chr == "Y") Chr <- 24
      if(!is.na(Chr) & Chr == "X, Y") Chr <- 23
      Chr <- as.numeric(Chr)
      RefSeq <- gsum[20]
      Chr.start <- as.numeric(gsum[21])
      Chr.end <- as.numeric(gsum[22])
      if(length(gsum)!=27){
        RefSeq <- gsum[19]
        Chr.start <- as.numeric(gsum[20])
        Chr.end <- as.numeric(gsum[21])
        }
      #Genom.start <- Chr.start
      #Genom.end <- Chr.end
    }
    gene.annot <- rbind(gene.annot, c(gname, Symb, Name, Org, Chr, Cytoband, Chr.start, Chr.end, RefSeq, Id))
  }	  
  gene.annot <- as.data.frame(gene.annot)
  colnames(gene.annot) <- c("Query", "Symb", "Name", "Org", "Chr", "Cytoband", "Chr.start", "Chr.end","RangeGB.Id", "GeneId")
  #	close(pb)
  cat('\n')
  return(gene.annot)
}

# gSearch: return gene Ids
gSearch <- function (geneSymb, database){  
  # geneSymb : use official symbols. Multiple requests are accepted, e.g. "EGFR, Homo sapiens"
  # database : have a look at 'http://eutils.ncbi.nlm.nih.gov/entrez/query/static/eutils_help.html' for details
  # 		on available databases and other e-tools as well.
  # ! This function can return more than one Id !
  
  gsrch.stem <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
  gsrch.mode <- paste("db=", database, "&retmode=xml","&term=", sep = "")
  URL <- paste(gsrch.stem, gsrch.mode, geneSymb, sep = "")
  doc <- xmlTreeParse(URL, isURL = TRUE, useInternalNodes = TRUE)
  sapply(c("//Id"), xpathApply, doc = doc, fun = xmlValue)
}

# gSummary: return gene summary
gSummary <- function (id, database){ 
  # id is provided by gsearch()
  # database : let's have a look at 'http://eutils.ncbi.nlm.nih.gov/entrez/query/static/eutils_help.html' for details
  # 		on available databases and other e-tools as well.
  
  sum.stem <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
  sum.mode <- paste("db=", database, "&id=", sep = "")
  doc <- xmlTreeParse(paste(sum.stem, sum.mode, id, sep = ""), isURL = TRUE, useInternalNodes = TRUE)
  sapply(c("//Item"), xpathApply, doc = doc, fun = xmlValue)
}

# http://eutils.ncbi.nlm.nih.gov/entrez/eutils/egquery.fcgi

# esearch: return pubmed Ids
eSearch <- function (term, n){
  # term = keywords. Same form as in a PubMed query
  # n = max number of pubmed Ids (papers) to return
  srch.stem <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
  srch.mode <- paste('db=pubmed&retmax=', n, '&retmode=xml&term=', sep = "")
  doc <-xmlTreeParse(paste(srch.stem,srch.mode,term,sep=""), isURL = TRUE, useInternalNodes = TRUE)
  sapply(c("//Id"), xpathApply, doc = doc, fun = xmlValue)
}

# myPubmed: build an html page containing papers : title, journal, date, abstract
myPubmed <- function(myKeywords, nAbst = 15, filename = NULL, title = NULL){
  pmids <- eSearch(myKeywords, nAbst)
  if(length(unlist(pmids)) == 0)
    stop ('No PubMed entries for:', myKeywords, '\n')
  x <- pubmed(pmids)
  a <- xmlRoot(x)
  numAbst <- length(xmlChildren(a))
  cat('kept', length(unlist(pmids)))
  arts <- vector("list", length = numAbst)
  for(i in 1:numAbst){
    arts[[i]] <- buildPubMedAbst(a[[i]])
  }
  if(is.null(filename)) filename = 'myPapers'
  filename <- paste(filename, '.html', sep = '')
  if(is.null(title)) title = myKeywords
  mypmHTML2(arts, filename = filename, title = title)
}

# Anchors: create links to the individual html pages
Anchors <- function(myTab, keyWords = NULL, nAbst= 15, inTitlAbst = TRUE){
  # create a new directory for web pages
  currentWd = getwd()
  htmlDocs = file.path(currentWd, 'htmlDocs')
  dir.create(htmlDocs)
  setwd(htmlDocs)

  # create html page with pubmed answers
  linkNames = tagNames = c()
  for(kw in myTab$Symb){
    cat('\n\n', kw, '\n')
    kwd1 = kw
    if(inTitlAbst) kwd1 <- paste(kw, "[Title/Abstract]", sep = '')
    myKeywords = paste(kwd1, "AND", keyWords)
    test = try(myPubmed(myKeywords, filename = paste(kw, keyWords), nAbst = nAbst), silent = TRUE)
    tmpName = NA
    if(is.null(test)) tmpName = paste(kw, keyWords)
    linkNames = c(linkNames, tmpName)
    tagNames = c(tagNames, kw)
  }

  # create links to html pages. These links will be added in the final table.
  anchors = ifelse(!is.na(linkNames),
                   paste("<a href=\'htmlDocs/", linkNames, '.html\'', " target=_blank >", tagNames, "</a>", sep=""),
                   'None')
  
  # back to the previous directory
  setwd(currentWd)
  return(anchors)
}

# Anchors: create links to the individual html pages
AnchorsSynapse <- function(parentId, myTab, keyWords = NULL, nAbst= 15, inTitlAbst = TRUE){

	# create html page with pubmed answers
  	linkNames = tagNames = c()
  	for(kw in myTab$Symb){
    	cat('\n\n', kw, '\n')
    	kwd1 = kw
    	if(inTitlAbst) kwd1 <- paste(kw, "[Title/Abstract]", sep = '')
    	myKeywords = paste(kwd1, "AND", keyWords)
    	test = try(myPubmed(myKeywords, filename = paste(kw, keyWords), nAbst = nAbst), silent = TRUE)
    	tmpName = NA
    	if(is.null(test)) tmpName = paste(kw, keyWords)
    	linkNames = c(linkNames, tmpName)
    	tagNames = c(tagNames, kw)
  		}		

  # create links to html pages. These links will be added in the final table.
  anchors = ifelse(!is.na(linkNames),
                   paste("<a href=\'htmlDocs/", linkNames, '.html\'', " target=_blank >", tagNames, "</a>", sep=""),
                   'None')
  # back to the previous directory
  setwd(currentWd)
  return(anchors)
}


# myHTML3: build the main html page
myHTML3 <- function(myTab, LogThres = 2, fileName = 'htmlTab', titleName = NULL){ 
  # Check for system
  # system <- Sys.info()["sysname"]
  # workingDir <- getwd()
  
  # Define links
  linkModel <- list()
  linkModel$to_Gene = "http://www.ncbi.nlm.nih.gov/sites/entrez?db=gene&cmd=Retrieve&dopt=Graphics&list_uids="					# Link to EntrezGene
  linkModel$to_Kegg = "http://www.genome.jp/dbget-bin/www_bget?hsa:"  															                      # Link to Kegg database
  
  # Define table title
  if(is.null(titleName)) titleName = fileName
  
  ## html	table initialization	
  html_init(init = T, filedir = getwd(), filename = fileName, titlename = titleName)	# changer titlename = tabName
  Append <- T
    
  ## Build complete links
  linkList <- linkModel
  linkList$to_Gene <- http_link(linkList$to_Gene, myTab$GeneId)
  linkList$to_Kegg <- http_link(linkList$to_Kegg, myTab$GeneId)
  
  # Formate table
  linkList2 <- c()
  for(L in 1:length(linkList)){
    tmplist <- sapply(as.data.frame(linkList[[L]]), function(x){paste("<a href=", x, " target=_blank >", myTab$GeneId, "</a>", sep="")})
    linkList2 <- cbind(linkList2, tmplist)
    }
  linkList2 <- as.data.frame(linkList2)
  colnames(linkList2) <-  names(linkModel)
  myTab = myTab[,-which(colnames(myTab) == 'GeneId')]
  if (!is.null(myTab$LogR)){
    Values = myTab$LogR
    styleL = function(x){paste("<p style=color:#0000FF;font-size:18;font-weight:bold;text-align:center>", x, "</p>", sep="")}
    styleN = function(x){paste("<p style=color:grey;font-size:18;text-align:center>", x, "</p>", sep="")}
    styleH = function(x){paste("<p style=color:#E62E00;font-size:18;font-weight:bold;text-align:center>", x, "</p>", sep="")}
    myTab$LogR  = ifelse(Values<(-LogThres), styleL(Values), ifelse(Values>LogThres, styleH(Values), styleN(Values)))
    }
  linkList2 <- cbind.data.frame(Symbol = myTab[,2], linkList2, myTab[,-c(1, 2)])  	# Information first(see above for #col concordance), then links.
  
  # Fill in html table
  html_tab(filename = fileName, data = linkList2)
  
  # Formate final html table 
  html_css(filename = fileName, css_liste = T)
  html_init(end = T, filedir = getwd(), filename = fileName)
  cat('file saved as: ', fileName, '.html in ', getwd(), sep = "")
}

# html.css: ccs values
html_css <- function(filename, css_liste){
  filename <-  paste(getwd(), '/', filename, ".html", sep = "")
  write.table("<style type=\"text/css\">", file = filename, quote=F, append=T, col.names=F, row.names=F)
  write.table("table{width:100%;}", file = filename, quote=F, append=T, col.names=F, row.names=F)
  write.table("tr.firstline{background-color:#FFBD9D;}", file = filename, quote=F, append=T, col.names=F, row.names=F)    		# fond entêtes de sous-tables
  write.table("a:link{text-decoration:none;color:blue;}", file = filename, quote=F, append=T, col.names=F, row.names=F)				# couleur du lien
  write.table("a:visited{text-decoration:none;color:#8A008A;}", file = filename, quote=F, append=T, col.names=F, row.names=F)			# couleur du lien après activation
  write.table("a:hover{text-decoration:underline;color:red;}", file = filename, quote=F, append=T, col.names=F, row.names=F)			# couleur du lien au passage de souris
  write.table("h2{background-color:#FFA366;text-align:center;}", file = filename, quote=F, append=T, col.names=F, row.names=F)		# fond entête principal
  write.table("span{font-weight:bold;}",file=filename, quote=F, append=T,col.names=F, row.names=F)
  write.table("</style> ",file = filename, quote=F, append=T,col.names=F, row.names=F)
}

#http_link
http_link <- function(http, id){
  http <- sapply(http, function(x){paste(x, id, sep = "")})
  return(http)
}

#html.init
## Initialize a html page.
# filedir : 	path
# filename: 	filename.
# titlename: 	main table title
html_init <- function(init = F, end = F, filedir = getwd(), filename = fileName, titlename = titleName){
  path <- paste(filedir, "/", filename, ".html", sep="")
  if(init){
    if(!file.exists(filedir)) system(paste("mkdir",filedir,sep=" "))
    HTMLInitFile(filedir, filename = filename, Title = titlename)
    HTML("<meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\" />",file=path)
    HTML(as.title(titlename),file=path)
    }
  if(end)	HTMLEndFile(file=path)
}

## Insert an image in a existing html page.
html_img <- function(filename, img_path){
  HTMLInsertGraph(GraphFileName=img_path, Width="50%", file=filename)
}

## Ajout de donnees contenues dans data Ã  une page html existante.
##filename: nom du fichier.	Argument obligatoires!
#data: donnÃ©es Ã  inclure dans le fichier
#caption: titre du tableau de donnees (valable si data est un dataframe/matrix)

#html_tab
html_tab <- function(filename, data=NULL, caption=NULL){
  path <- paste(filename, ".html", sep = "")
  if(!is.null(data)){HTML(data, row.names = F, innerBorder = 1, caption = caption, captionalign = "top", file = path)}
}

#mypmHTML2
mypmHTML2 <- function (absts, filename, title, table.center = TRUE) {
  frames = FALSE
  if (!is.list(absts)) {
    if (is(absts, "pubMedAbst")) 
      absts <- list(absts)
    else stop("'absts' parameter does not seem to be valid.")
  }
  if (missing(filename))
    filename <- "absts.html"
  if (missing(title)) 
    title <- "BioConductor Abstract List"
  nrows = length(absts)
  pmids <- unlist(lapply(absts, pmid))
  dates <- unlist(lapply(absts, pubDate))
  journals <- unlist(lapply(absts, journal))
  abstrs <- unlist(lapply(absts, abstText))
  queries <- unlist(lapply(absts, 
                           function(x) {
                             pm <- pmid(x)
                             out <- pmidQuery(pm)
                             out
                           })
                    )
  titles <- unlist(lapply(absts, articleTitle))
  # anchors <- makeAnchor(queries, titles, toMain = frames)
  anchors <- makeAnchor2(queries, titles)
  topText <- paste("<html>\n<head>\n<title>", title, "</title>", 
                   "\n</head>\n<body bgcolor=#E6E2AF>\n", "<H1 ALIGN=CENTER>", 
                   title, "</H1>\n", "</body></title>", sep = "")
  
  head <- c("Article Title", "Journal", "Publication Date", "Abstract")
  headOut <- paste("<TH>", head, "</TH>", collapse = "\n")
  
  outfile <- file(filename, "w")
  cat(topText, file = outfile)
  if (table.center) 
    cat("<CENTER> \n", file = outfile)
  cat("<TABLE BORDER=1>", file = outfile, sep = "\n")
  cat("<TR bgcolor=#FFB894>", headOut, "</TR>", file = outfile, sep = "\n")
  tds <- paste("<TD width=445, bgcolor=#F5F3DF>", anchors, 
               "</TD><TD width=200, bgcolor=#F5F3DF>", journals, 
               "</TD><TD width=135, bgcolor=#F5F3DF>", dates,
               "</TD><TD bgcolor=#F5F3DF>", abstrs, "</TD>",
               sep = "", collapse = "\n</TR>\n<TR>\n")
  for (td in tds) cat("<TR>", td, "</TR>", file = outfile, sep = "\n")
  cat("</TABLE>", file = outfile)
  if (table.center) 
    cat("</CENTER> \n", file = outfile)
  cat("</body>", "</html>", sep = "\n", file = outfile)
  close(outfile)
  invisible(NULL)
}

# Replace annotate::makeAnchor
makeAnchor2 <- function (link, title) {
  Target = '_blank>'
  out <- paste('<A HREF=', link, ' target=', Target, title, '</A>', sep = "")
  out
}

