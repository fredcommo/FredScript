# Build RDSobject from CEL files and samples data

require(affy)
require(synapseClient)
# args <- commandArgs(True)
# parentId <- args[1]
# CELid <- args[2]
# sampleId <- args[3]

buildRDS <- function(parentId){
  cat('building', parentId, '\n')
  query <- synapseQuery(paste0('select name, id from entity where parentId =="', parentId, '"'))
  CELid <- query$entity.id[grep('CELfiles', query$entity.name)]
  sampleId <- query$entity.id[grep('[Ss]amples', query$entity.name)]
  
  cat('Loading data from synapse...\t')
  CEL <- synGet(CELid)
  entSamples <- synGet(sampleId)
  cat('Done\n')
  
  cat('Reading samples...\n')
  samples <- read.csv(entSamples@filePath, header = TRUE, sep = '\t')
  renameCol <- lapply(1:ncol(samples), function(i){
    if(any(grepl(':', samples[,i]))){
      values <- as.data.frame(gsub('^[a-zA-Z ]+\\: ', '', samples[,i]))
      item <- as.character(samples[1,i])
      r <- regexpr('^[a-zA-Z ]+', item)
      colName <- regmatches(item, r)
      colnames(values) <- gsub(' ', '_',colName)
    }
    else{
      values <- as.data.frame(samples[,i])
      colnames(values) <- colnames(samples)[i]
    }
    cat(i, colnames(values), '\n')
    return(values)
  })
  samples <- do.call(cbind.data.frame, renameCol)
  samples <- samples[order(samples$Sample_geo_accession),]
  cat('Done\n')
  
  cat('Reading CEL files...\t')
  celPath <- paste0(gsub('NOT_SET', '', CEL@filePath), 'CEL')
  dir.create(celPath, recursive = TRUE)
  untar(CEL@filePath, exdir = celPath)
  setwd(celPath)
  rawData <- ReadAffy()
  cat('Done\n')
  
  cat('Normalize...\n')
  eset <- rma(rawData)
  eset <- exprs(eset)
  eset <- eset[ ,order(colnames(eset))]
  arrayIds <- gsub('\\.(.*)+', '', colnames(eset))
  sampleIds <- samples$Sample_geo_accession
  if(ncol(eset)!=nrow(samples))
    cat('Attention: number of arrays do not correspond to number of samples\n' )
  if(!all(arrayIds == sampleIds))
    cat('Attention: arrays and samples do not match\n' )
  cat('Done\n')
  
  cat('Computing A/P calls...\n')
  esetCall <- exprs(mas5calls(rawData))
  freqCall <- lapply(1:nrow(esetCall), function(i){
    A <- sum(esetCall[i,]=='A')
    M <- sum(esetCall[i,]=='M')
    P <- sum(esetCall[i,]=='P')
    c(A = A/ncol(esetCall), M = M/ncol(esetCall), P = P/ncol(esetCall))
  })
  freqCall <- do.call(rbind, freqCall)
  freqCall <- apply(freqCall, 1, function(x){names(which.max(x))})
  freqCall <- as.data.frame(freqCall)
  rownames(freqCall) <- rownames(esetCall)
  cat('Done\n')
  
  #########################
  # Store in synapse
  cat('Pushing to synapse...\n')
  fileName <- paste0(parentId, '.rds')
  Data <- list(eset = eset, samples = samples, freqCall = freqCall)
  save(Data, file = file.path(tempdir(), fileName))
  file <- File(file.path(tempdir(), fileName), parentId = parentId)
  file <- synStore(file,
                   activityName = paste('CEL files preprocessing -', rawData@annotation, '- Human samples'),
                   used = list(list(entity = CEL, wasExecuted = FALSE),
                               list(entity = entSamples, wasExecuted = FALSE)))
  # Update the wiki
  cat('Updating the wiki...\n')
  folderWikiUri <- sprintf("/entity/%s/wiki", parentId)
  folderWiki <- synRestGET(folderWikiUri)
  folderUpdateWikiUri <- sprintf("/entity/%s/wiki/%s", parentId, folderWiki$id)
  folderWiki$attachmentFileHandleIds[length(folderWiki$attachmentFileHandleIds)+1] <- file@fileHandle$id
  newWiki <- paste0('###', parentId,
                    '.rds contains both preprocessed values (Log2) and samples information\n',
                    'Version1 is normalized using mas5\nVersion2 is normalized using rma\n\n')
  folderWiki$markdown <- paste0(newWiki, '\n', folderWiki$markdown,'\n')
  folderWiki <- synRestPUT(folderUpdateWikiUri, folderWiki)
  #  unlink(synapseCacheDir(), recursive = TRUE, force = TRUE)
  #  unlink(tempdir(), recursive = TRUE, force = TRUE)
  cat('Done\n')
  #########################
}

#'syn1935130', 'syn1939680', 'syn1935134', 'syn1939658', 'syn1939654', 'syn1939656'
synIds <- c('syn1939674', 'syn1939669', 'syn1939676', 'syn1939671', 'syn1939678')
buildRDS('syn1939656')
