
# The function activates an interactive selection of file to load
# Collect microarray informations and Cy5/Cy3 intensity values,
# and suppress flags and duplicated probes.
# Root : A letter indicating the hard drive partition to use. May run with other indications (~/home)

##########################################################################################
###########################

# Class are defined in AllClasses.R
# Accessors 	are defined in AllAccessors.R
# Generic methids are defined in AllGenerics.R
# Methods are in AllMethods .R

# validity method : To do !

# Common functions : To do !
	# setGeneric('Norm', function(object) standardGeneric('Norm'))
	# setGeneric('Segment', function(object) standardGeneric('Segment'))
	# setGeneric('Plot', function(object) standardGeneric('Plot'))
	# setGeneric('Save', function(object) standardGeneric('Save'))
	# setGeneric('EditSupTab', function(object) standardGeneric('EditSubTab'))



###########################
###########################

buildCGHObj.v01 <- function(){
	f = try(getFile(), silent = TRUE)								# helper function
	if(class(f) == 'try-error') stop ('No selected file!\n')
	
	cat("\nGetting Array Information...")
	object = getAnnot(f)												# helper function
	object <- readInfo(object)										# Class function
	
	cat('\n\nReading', f, '...')
	object <- readCN(object)											# Class function
	
	cat('\n\nCleaning data...')	
	object <- supressFlags(object)									# Class function
	object <- supressDuplic(object)								# Class function
	object <- preset(object)											# Class function
	
	cat("\n\n\t*** Everything is fine. Lucky you :-) ***\n\n")
	return(object)
}


###########################
getFile <- function(){
	
	'
	Called by buildCGHObj()
	Open a dialogue box and return the file selected by the user.
	'

	#fileName <- tclvalue(tkgetOpenFile()) 									# Open current the directory to select the file to load.
	fileName = file.choose()
	if (!nchar(fileName)) 
		tkmessageBox(message = "No file selected!")
	#else 
    #	tkmessageBox(message = paste("Selected File:", fileName))

	fileName = unlist(strsplit(fileName, '/'))
	path = fileName[1]
	for (i in 2:(length(fileName)-1))
		path = paste(path, fileName[i], sep = '/')
	setwd(path)
	return (fileName[length(fileName)])															# add Folder
}

###########################
getAnnot <- function(fileName){
	'
	Called by buildCGHObj()
	Open the selected file and read the first line to identify the type of array (platform = Agilent or Affy)
	Return an object of class Agilent or Affymetrix.
	Assign the fileName (without its path), the sampleId, and the platform to object@info
	'	
	obj = list()
	fileInfo = unlist(strsplit(fileName, '_'))
	sampleId = 	fileInfo[1]
	a <- try(read.csv(fileName, header = F, fill = T, skip = 0, nrows = 1, sep = "\t"), silent = T)
	{
		if(class(a)!="try-error") cat("\nFile summary:", "\n\tpath: ", getwd(), "\n\tfileName:", fileName)
		else stop("\t", fileName, ": No such file on directory (glup!). Please fix it!\n\n")
	}

	# According to 'platform', create a corresponding cghObject
	if (a[1,1] == '#GenomeWideSNP_6.na32.annot.db'){
		platform = 'Affymetrix'
		Obj <- AffyObj(info = c(fileName = fileName, sampleId = sampleId, platform = platform))
		}
	else{
		platform = 'Agilent'
		Obj <- AgilentObj(info = c(fileName = fileName, sampleId = sampleId, platform = platform))
	}
	return (Obj)
}

###########################

CyAdjust <- function(cnSet, Fract, Centr){
	'
	Called by adjustSignal(object)
	If Fract (prop of tumor in the sample), adjust the tumor signal (Cy5)
	Compute and adjust the Log2(Cy3/Cy5)
	Return the data.frame with a supplementary column: Log2Ratio
	'
	cat('\nCy effect adjustment...')
	g <- log2(cnSet$gMedianSignal)					# Ref in Cy3
	r <- log2(cnSet$rMedianSignal)					# Test in Cy5

	# Calculates weights to correct the dilution effect (due to tumor cell rate). No effect if Fract = 1. DO NOT USE until validation !
	if(!is.null(Fract)){
		Q <- quantile(r, Fract, na.rm = TRUE)
		w <- 1/(1+exp(1/sqrt(Fract)*(Q-r)))
		r <- r*(1 + (1-Fract)*w)
		}
	M <- r - g
	A <- (r + g)/2
	Loess <- loessFit(M, A)$fitted
	LR <- M - Loess
	cnSet$Log2Ratio <- LR 	
	if (Centr)
		cnSet$Log2Ratio - median(cnSet$Log2Ratio, na.rm = T)
	cat('\tDone.')
	return (cnSet)
	}

#################

GCadjust <- function(cnSet){
	'
	Called by adjustSignal(object)
	Adjust the GC%
	'
	arrayInfoPath = '/gluster/home/jcommo/ArraysInfos/'
	cat('\nGC% adjustment...')
	AgilentDB = readRDS(paste0(arrayInfoPath, 'Agilent_022060_4x180_hg19_20130116_FC.RData'))
	AgilentDB <- AgilentDB[which(AgilentDB$ProbeID %in% cnSet$ProbeName),]
	# Check the probeNames
	if(!all(as.character(AgilentDB$ProbeID) == as.character(cnSet$ProbeName)))
		stop('Agilent DB: Probe names do not match.')
	lr = cnSet$Log2Ratio
	GC <- AgilentDB$GCpercent
	adjLr <- lr - loessFit(lr, GC)$fitted
	cnSet$Log2Ratio = adjLr
	cat('\tDone.')
	cnSet <- cbind.data.frame(cnSet[,c('ProbeName', 'ChrNum', 'ChrStart')], genomicPos = AgilentDB$genomicPos, Log2Ratio = cnSet[,'Log2Ratio'])
	return(cnSet)
	}

#################

smoothLR <- function(LR, Platform, cut, K){
	'
	Called by EMnormalize(object)
	Smoothing the LR vector to improve the EM classification.
	'
	if(Platform == 'Affymetrix')
		LR <- LR[seq(1, length(LR), by = 6)]
	if(any(is.na(LR))) LR <- LR[!is.na(LR)]
	runLR <- runmed(LR, k = K)	
	q1 <- cut[1]; q2 <- cut[2]
	runLR = runLR[which(runLR>=q1 & runLR<=q2)]
	return (runLR)
}

#################

buildEMmodel <- function(LR, G, by){
	'
	Called by EMnormalize(object)
	Model the distribution as a gaussian mixture.
	'
	model <- Mclust(LR[seq(1, length(LR), by = by)], G = G)
	nG <- model$G
	p <- model$parameters$pro
	m <- model$parameters$mean
	s <- model$parameters$variance$sigmasq
	if(length(s)<length(m)) s <- rep(s, length(m))
	p <- p[order(m)]
	s <- s[order(m)]
	m <- m[order(m)]
	return(list(nG = nG, m = m, p = p, s =s))
}

#################

mergePeaks <- function(nG, m, s, p, MergeVal){
	'
	Called by EMcentr(object)
	Given the initial EM parameters, recompute the gaussian mixture by merging the groups for which the between-group distance < MergeVal
	'
	Raw <- c(1, 1)
	while(length(Raw)!=0){
		Mdist <- matrix(0, nG, nG)
		for(i in 1:nG)
			for(j in 1:nG){
				Mdist[i, j] <- abs(m[i] - m[j])
				}
			diag(Mdist) <- NA
			Raw <- ceiling(which(Mdist<MergeVal)/nG)
			cat('Merging', Raw, "\n")
			if(length(Raw)!=0){
				C1 <- Raw[1]
				C2 <- Raw[2]
				m[C1] <- ((p[C1])*m[C1] + (p[C2])*m[C2])/(p[C1] + p[C2])
				s[C1] <- ((p[C1])*s[C1] + (p[C1])*s[C2])/(p[C1] + p[C2])
				p[C1] <- p[C1] + p[C2]
				m <- m[-C2]; s <- s[-C2]; p <- p[-C2]
				nG <- length(m)
				cat("means:", m, "\nVar:", s, "\nprops:", p, "\n\n")
				}
			}
		return(list(nG = nG, m = m, s = s, p = p))
		}


#################

computeDensities <- function(n, m, p, s){
	'
	Called in EMnormalize(object)
	Simulates the mixture model according to the returned EM paramaters.
	'
	dList = list()
	peaks <- kurt <- c()
	for(i in 1:length(m)){
		tmp <- rnorm(n*p[i], m[i], sqrt(s[i]))
		tmpD <- density(tmp, na.rm = T)
		tmpD$y = tmpD$y *p[i]
		dList[[i]] <- tmpD
		kurt <- c(kurt, kurtosis(tmp, type = 2))
		peaks <- c(peaks, max(tmpD$y))
		}
	return(list(dList = dList, peaks = peaks, kurt = kurt))
}

#################

chooseBestPeak <- function(peaks, m, peakThresh){
	'
	Called in EMnormalize(object)
	Estimates what peak as to be used as the centralization value.
	'
	best <- which(peaks>=max(peaks)*peakThresh & m<=0)
	if (length(best) != 0)
		cat('\nLeft peak at', m[best], 'has been chosen.')
	else{
		# Centered peak next
		best <- which(peaks>=max(peaks)*peakThresh)
		best <- best[which.min(abs(m)[best])]
		if (length(best) != 0)
			cat('\nCentral peak at', m[best], 'has been chosen.')
		else{
			# Right peak last
			best <- which(peaks>=max(peaks)*peakThresh & m>=m[which.max(peaks)])
			if (length(best) != 0)
					cat('\nRight peak at', m[best], 'has been chosen.\n')
				}
			}
	bestPeak = best[which.max(m[best])]
	return(bestPeak)
}

#################

plotEMmodel <- function(LR, dList, m, bestPeak, Title){
	'
	Called in EMnormalize(object)
	Visualization of the mixture model
	'
	dLR <- density(LR)
	currentPlot = xyplot(dLR$y~dLR$x, type = "n", main = Title,
						xlab = "Log2R", ylab = "Density",
						xlim = range(-0.6, 0.6), ylim = range(0, max(dLR$y)*1.5),
						panel = function(x, y){
											lpolygon(dLR$x, dLR$y, col = 'grey90')
											n = length(LR)
											nG = length(dList)
											for (i in 1:nG){
												tmp <- dList[[i]]
												llines(tmp$x, tmp$y, lwd = 1, col = rgb(i/nG, 0.2, (nG-i)/nG, 0.75))
												lpolygon(tmp$x, tmp$y, col = rgb(i/nG, 0.2, (nG-i)/nG, 0.25))
												ltext(	x = mean(tmp$x), y = min(max(tmp$y*1.5, na.rm = TRUE), max(dLR$y)*1.25), labels = round(m[i], 3),
														cex = ifelse(i == bestPeak, 1.5, 1.25), font =  ifelse(i == bestPeak, 2, 1))
												}
									}
					)
	return(currentPlot)
}

# To Do

#################
# QCSegm

#################
MedSegm <- function(seg.cna.obj, cnSet){
	'
	Called by SegmentCGH()
	'
	seg.start <- seg.cna.obj$output$loc.start
	seg.end <- seg.cna.obj$output$loc.end			
	seg.len <- seg.cna.obj$output$num.mark
	cum.seg.len <- cumsum(seg.len)
	s <- cnSet$ChrStart[cum.seg.len]
	e <- cnSet$ChrEnd[cum.seg.len]
	if(!is.null(e)) seg.end <- seg.end + (e-s)


	seg.med <- c()
	for(i in 1:length(seg.start)){
		index <- which(cnSet$genomicPos>= seg.start[i] & cnSet$genomicPos<= seg.end[i])
		tmpLR = cnSet$Log2Ratio[index]
		tmpMed <- tukey.biweight(tmpLR[!is.na(tmpLR)])
		seg.med <- c(seg.med, tmpMed)
		}
	seg.cna.obj$output <- cbind.data.frame(seg.cna.obj$output, seg.med = seg.med)
	return(seg.cna.obj)
}

#################
AddSegments <- function(seg.cna.obj, cnSet, cutGainLoss = 10, use.medians = TRUE){
	'
	Called by SegmentCGH()
	'
	seg.start <- seg.cna.obj$output$loc.start
	seg.end <- seg.cna.obj$output$loc.end			# !!!! verifier seg-end
	seg.mean <- seg.cna.obj$output$seg.mean
	if(use.medians) seg.mean <- seg.cna.obj$output$seg.med
	seg.len <- seg.cna.obj$output$num.mark
	cum.seg.len <- cumsum(seg.len)
	s <- cnSet$ChrStart[cum.seg.len]
	e <- cnSet$ChrEnd[cum.seg.len]
	if(!is.null(e)) seg.end <- seg.end + (e-s)

	# Proportion d'aberrations
	cutoff <- log2(1 + cutGainLoss/100)
	seg.delta <- seg.end - seg.start
	Prop <- sum(seg.delta[which(abs(seg.mean)>=cutoff)])/sum(seg.delta)*100
	Prop <- round(Prop, 2)

	seg.val <- rep(NA, nrow(cnSet))
	for(i in 1:length(seg.mean)){
		index <- which(cnSet$genomicPos>= seg.start[i] & cnSet$genomicPos<=seg.end[i])
		seg.val[index] <- seg.mean[i]
		}
	cnSet <- cbind.data.frame(cnSet, Segm = seg.val)
	return(cnSet)
	}

#################
dLRsd <- function(LR){
	'
	Not used yet
	'
	n <- length(LR)
	V1 <- LR[-1]
	V2 <- LR[-n]
	dLR <- V2-V1
	q1 <- quantile(dLR, 0.25, na.rm = TRUE)
	q3 <- quantile(dLR, 0.75, na.rm = TRUE)
	s <- sd(dLR[which(dLR > q1 & dLR < q3)], na.rm = TRUE)/sqrt(2)
	return(s)
	}

#################
# Gene list functions
#################

# http://eutils.ncbi.nlm.nih.gov/entrez/eutils/egquery.fcgi

# esearch: return pubmed Ids
eSearch <- function (term, n){
	'
	Not used yet
	'
  # term = keywords. Same form as in a PubMed query
  # n = max number of pubmed Ids (papers) to return
  srch.stem <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
  srch.mode <- paste('db=pubmed&retmax=', n, '&retmode=xml&term=', sep = "")
  doc <-xmlTreeParse(paste(srch.stem,srch.mode,term,sep=""), isURL = TRUE, useInternalNodes = TRUE)
  sapply(c("//Id"), xpathApply, doc = doc, fun = xmlValue)
}

# gSearch: return gene Ids
gSearch <- function (geneSymb, database){  
  '
	Called by geneRequest.v7()
  '
	require(XML)
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

itemList = c(	'Orgname', 'NomenclatureStatus', 
              'NomenclatureSymbol', 'OtherAliases', 'Description',
              'Chromosome', 'MapLocation', 'ChrStart', 'ChrStop')

gSummary <- function(id, database, itemList){
  '
	Called by geneRequest.v7()
  '
	require(XML)
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


geneRequest.v7 <- function(geneId, database, verbose = TRUE){
  	output = data.frame()
  	notFound = cbind.data.frame(NA, t(rep(NA, length(itemList))))
    colnames(notFound) = c(itemList, 'entrezgene')
    if(is.character(geneId)){
    	geneId = toupper(geneId)
	    id = gSearch(paste0(geneId, '[symbol]+homo+sapiens[Organism]'), database)
	    id = unlist(id)
	    }
    if(length(id) == 0){
      if(verbose) cat('\n', geneId, '\t*** not found ***')
      output = rbind.data.frame(output, notFound)
    	}
    else{
      # should have NomenclatureStatus = 'Official'
      Official = FALSE
      k = 1
      while (!Official & k <= length(id)){
        tmp = gSummary(paste0(id[k], '%5BUID%5D'), database, itemList)
        tmp = cbind.data.frame(tmp, entrezgene = id[k])
        Official <- tmp$NomenclatureStatus == 'Official'
        k = k + 1
      	}
      if(!Official) tmp = notFound
      output = rbind.data.frame(output, tmp)
      if(verbose){
      	if (is.character(geneId)) cat('\n', geneId, 'found. entrezgene:', as.character(output$entrezgene))
      	else cat('\n', geneId, 'found as', as.character(output$NomenclatureSymbol))
      	}
    }
 
	# Add genomic position.
	arrayInfoPath = '/gluster/home/jcommo/ArraysInfos/'
	hg19 <- read.csv(paste(arrayInfoPath, 'human.chrom.info.hg19.FC.txt', sep = ''), header = TRUE, sep = '\t')
	cumLen = cumsum(as.numeric(hg19$length))
	cumLen = c(0, cumLen[-length(cumLen)])
	output = cbind.data.frame(output[,-ncol(output)], genomicStart = rep(NA, nrow(output)), genomicStop = rep(NA, nrow(output)), 														entrezgene = output[,ncol(output)])
	chr = as.numeric(as.character(output$Chromosome))
	output$genomicStart = as.numeric(as.character(output$ChrStart))+ cumLen[chr]
	output$genomicStop = as.numeric(as.character(output$ChrStop)) + cumLen[chr]
  	#cat('\n')
  	#rownames(output) = seq(1, nrow(output))
  	return(output)
}


# geneRequest <- function(geneList){
	# require(biomaRt)
	# human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
	# output = getBM(attributes=c('hgnc_symbol', 'description', 'chromosome_name', 'band', 'start_position','end_position', 'entrezgene'),
			# filters = 'hgnc_symbol', values = geneList, mart = human)
	# output$description =  gsub(' \\[.*\\]', '', output$description)
	# output = as.data.frame(output)
	# # Add genomic positions
	# arrayInfoPath = '/Users/fredcommo/Documents/Projet Safir/Arrays Infos/'
	# hg19 <- read.csv(paste(arrayInfoPath, 'human.chrom.info.hg19.FC.txt', sep = ''), header = TRUE, sep = '\t')
	# cumLen = cumsum(as.numeric(hg19$length))
	# cumLen = c(0, cumLen[-length(cumLen)])
	# output = cbind.data.frame(output[,-ncol(output)], genomicStart = rep(NA, nrow(output)), genomicEnd = rep(NA, nrow(output)), entrezgene = output[,ncol(output)])
	# foreach(i = 1:nrow(output)) %do%{
		# chr = as.numeric(output$chromosome_name[i])
		# output$genomicStart[i] = output$start_position[i] + cumLen[chr]
		# output$genomicEnd[i] = output$end_position[i] + cumLen[chr]
		# }
	# return(output)
# }


#################
# Plot functions
#################

locateChr <- function(y){
	colText = 'grey40'
	colLines = 'grey80'
	arrayInfoPath = '/gluster/home/jcommo/ArraysInfos/'
	hg19 <- read.csv(paste(arrayInfoPath, 'human.chrom.info.hg19.FC.txt', sep = ''), header = TRUE, sep = '\t')
	cumLen = cumsum(as.numeric(hg19$length))
	cumCentr <- 1/2*cumLen[1]
	for(chr in 2:length(cumLen)) cumCentr = c(cumCentr, cumLen[chr-1] + 1/2*cumLen[chr])
	panel.abline(h = 0)#, lty = 3)
	panel.abline(v = cumLen[1:23], col = colLines, lty = 2)
	ltext(0, y, labels = "chr", cex = 0.75, col = colText)
	ltext(cumLen[1]/2, y, labels = 1, cex = 0.75, col = colText)
	for(i in 2:length(cumCentr)){
		x <- (hg19$length[i]/2 + cumLen[i-1])
		ltext(x, y, labels = i, cex = 0.75, col = colText)
		}
	}

addRunMed <- function(genomicPos, LogRatio){
	if(any(is.na(LogRatio))) LogRatio = LogRatio[!is.na(LogRatio)]
	Samp = seq(1, length(LogRatio), len = length(LogRatio)/10)
	g = genomicPos[Samp]; rmed = runmed(LogRatio[Samp], k = length(LogRatio)/1500)
	llines(g, rmed, col = 'black', lwd = 0.5)
	}

locateProbes <- function(genomicPos, LogRatio){
	if(any(is.na(LogRatio))) LogRatio = LogRatio[!is.na(LogRatio)]
	Samp = seq(1, length(LogRatio), len = length(LogRatio)/20)
	g = genomicPos[Samp]; LR = LogRatio[Samp]
	lpoints(g, LR, pch = 19, cex = 0.1, col = 'grey80')
	#rgb(0.8, 0.8, 0.8, 0.5)
	}

locateSegments <- function(genomicPos, Segments, thresh){
	Samp = seq(1, length(Segments), len = length(Segments)/10)
	g = genomicPos[Samp]; s = Segments[Samp]
	lpoints(g, s, pch = 19, cex = 0.25, col = ifelse(s>thresh, 'dodgerblue3', ifelse(s < -thresh , 'red3', 'grey50')))
	}


####################
#	Table functions
####################

keggLink <- function(id){
	return(paste0('http://www.genome.jp/dbget-bin/www_bget?hsa:', id))
}
ncbiLink <- function(id){
	return(paste0('http://www.ncbi.nlm.nih.gov/gene?term=', id, '%5BUID%5D'))
}
htmlLink <- function(tag, link){
	return(paste0("<font size='4'><a href=", link, " target=_blank >", as.character(tag), "</a></font>"))
}

html_css <- function(filename){
		write.table("<style type=\"text/css\">", file = filename, quote=F, append=T, col.names=F, row.names=F)
		write.table("table{width:100%;height:10px;background-color:#FBFBEF;}", file = filename, quote=F, append=T, col.names=F, row.names=F)
		write.table("tr.firstline{background-color:#FFBD9D;}", file = filename, quote=F, append=T, col.names=F, row.names=F)				# fond entêtes de sous-tables
		write.table("a:link{text-decoration:none;color:blue;}", file = filename, quote=F, append=T, col.names=F, row.names=F)				# couleur du lien
		write.table("a:visited{text-decoration:none;color:#8A008A;}", file = filename, quote=F, append=T, col.names=F, row.names=F)			# couleur du lien après activation
		write.table("a:hover{text-decoration:underline;color:red;}", file = filename, quote=F, append=T, col.names=F, row.names=F)			# couleur du lien au passage de souris
		write.table("h2{background-color:#FFA366;text-align:center;}", file = filename, quote=F, append=T, col.names=F, row.names=F)		# fond entête principal
		write.table("span{font-weight:bold;}",file=filename, quote=F, append=T,col.names=F, row.names=F)
		write.table("#Norm{color:#A1A0A0;}",file=filename, quote=F, append=T,col.names=F, row.names=F)
		write.table("#Gain{color:#3075ED;}",file=filename, quote=F, append=T,col.names=F, row.names=F)
		write.table("#Ampli{color:#1100A6;}",file=filename, quote=F, append=T,col.names=F, row.names=F)	
		write.table("#Loss{color:#E50810;}",file=filename, quote=F, append=T,col.names=F, row.names=F)												# 
		write.table("</style> ",file = filename, quote=F, append=T,col.names=F, row.names=F)
	}
	
textStyle <- function(x, align = 'center'){
	#return(paste0('<font size="4", align="center">', x, "</font>"))
	return(paste0("<p style=font-size:18;text-align:", align, ">", x, "</p>"))
}

logStyle<- function(LR, thresh){
	styleL = function(x){paste0("<p style=color:#E62E00;font-size:18;font-weight:bold;text-align:center>", round(x, 3), "</p>")}
    styleN = function(x){paste0("<p style=color:grey;font-size:18;text-align:center>", round(x, 3), "</p>")}
    styleH = function(x){paste0("<p style=color:#0000FF;font-size:18;font-weight:bold;text-align:center>", round(x, 3), "</p>")}
    LR  = ifelse(LR<(-thresh), styleL(LR), ifelse(LR>thresh, styleH(LR), styleN(LR)))
	return(LR)
}

buildHtml <- function(geneTable, filePath, fileName, cssFile){

	geneTable$Chromosome = as.numeric(as.character(geneTable$Chromosome))
	geneTable$ChrStart = as.numeric(as.character(geneTable$ChrStart))
	geneTable$ChrStop = as.numeric(as.character(geneTable$ChrStop))
	geneTable = geneTable[order(geneTable$Chromosome, geneTable$ChrStart),]
	isSymbol = grep('[Ss]ymbol', colnames(geneTable))
	desc = grep('[Dd]escription', colnames(geneTable))
	chr = grep('[Cc]hromosome', colnames(geneTable))
	mapLoc = grep('[Mm]ap[Ll]oc', colnames(geneTable))
	cStart = grep('[^genomic][Ss]tart', colnames(geneTable))
	cStop = grep('[^genomic][Ss]top', colnames(geneTable))
	toNcbi <- htmlLink(geneTable[,isSymbol], ncbiLink(geneTable$entrezgene))
	toKegg <- htmlLink(geneTable[,isSymbol], keggLink(geneTable$entrezgene))
	TitleName <- paste0("<h1 align=center>Genes of Interest</h1><h2 align=center>Sample ", fileName, '</h2>')
	output = HTMLInitFile(filePath, filename = fileName, Title = 'Gene Of Interest')
	geneTable = cbind.data.frame(Symbol = toNcbi,
												Description = textStyle(geneTable[,desc], align = 'left'),
												Chr = textStyle(geneTable[,chr]),
												MapLoc = textStyle(as.character(geneTable[,mapLoc])),
												ChrStart = textStyle(round(geneTable[,cStart]/1e3)),
												ChrEnd = textStyle(round(geneTable[,cStop]/1e3)),
												geneId = textStyle(geneTable$entrezgene),
												Log2Ratio = logStyle(geneTable$Log2Ratio, 0.1))
	colnames(geneTable)[5:6] = c('ChrStart (Kb)', 'ChrEnd (Kb)')
	HTML("<meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\" />", file = output)
	HTML(as.title(TitleName), file=output)
	HTML(geneTable, file = output, row.names = F, innerBorder = 1, caption = NULL, captionalign = "top", CSSFile = cssFile)
	html_css(output)
	HTMLEndFile()		
}


FullHtml <- function(segTable, cutoff){
	# use biomaRt here to get all the genes included between 2 positions.
	arrayInfoPath = '/Users/fredcommo/Documents/Projet Safir/Arrays Infos/'
	myBiomart <- readRDS(paste0(arrayInfoPath,'myBiomaRt.rds'))

}


