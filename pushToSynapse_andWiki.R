# demo of how to create a file then use the uploaded file in a wiki
projectId <- "syn1834037"

# note, we recommend using "File" rather than "Data"
myPath <- '/Users/fredcommo/Documents/MyProjects/Projet ACP/Some_TCGA/'
myImage <- "unc.edu_BLCA_IlluminaHiSeq_RNASeqV2.geneExp.whitelist_tumor_syn1571504_PCAfilt_2.png"

file <- File(paste0(myPath, myImage), parentId=projectId)
file <- synStore(file)
fileHandleId <- file@fileHandle$id
fileWikiUri <- sprintf("/entity/%s/wiki", propertyValue(file, "id"))

# we start a wiki
fileWiki<-list()

# we add to our wiki the ID of the previously uploaded file
fileWiki$attachmentFileHandleIds<-list(fileHandleId)

# in the markdown we say to display the image.  Note, 'fileName' is the URLEncoded version of the file chosen above.
fileWiki$markdown<-"${image?fileName=mygraph%2Epng}"

# now 'push' the wiki to Synapse
fileWiki<-synRestPOST(fileWikiUri, fileWiki)

# voila!
onWeb(propertyValue(file, "id"))


######################################
# Creates folders
# push images in each folder
# updates the folder's wiki page to visualize the images.
# FC/2013_05_08
######################################

#####################################
# Push to synapse
setwd('/Users/fredcommo/Documents/MyProjects/Projet ACP/Some_TCGA/')

Path <- paste0(getwd(), '/')
projectId = 'syn1834953'
myCode <- Code(list(name = "exploreTCGA_PCAfiltering.R", parentId = projectId))
myCode <- addFile(myCode, paste0(Path, 'exploreTCGA_PCAfiltering.R'))
myCode <- storeEntity(myCode)
synIds <- rep(rnaSeqList$entity.id, each = 2)


#####################################
# Push files & create a wiki page
#myCode <- getEntity('syn1834779')

listFiles = list.files()
for (l in 2:length(listFiles)){
  fileName = listFiles[l]
  #  synId = synIds[l-1]
  #  synId = rnaSeqList$entity.id[i]
  splitFileName <- unlist(strsplit(fileName, '_'))
  tissue <- splitFileName[2]
  synId <- splitFileName[6]
  cat(tissue, synId, '\n')
  
  folderName <- paste(tissue, synId, sep = '_')
  Folders <- synapseQuery('select id, name from entity where entity.parentId == "syn1834953"')
  if(!folderName %in% Folders$entity.name){
    # Create a new folder
    newFolder <- Folder(name = folderName, parentId = projectId)
    newFolder <- synStore(newFolder)
    parentId <- propertyValue(newFolder, 'id')
    # Create the folder wiki uri.    
    folderWikiUri <- sprintf("/entity/%s/wiki", parentId)
    # Start a wiki, and initialize the markdown as empty.
    folderWiki <- list()
    folderWiki$attachmentFileHandleIds <- list()
    folderWiki$markdown <- c()
    folderWiki <- synRestPOST(folderWikiUri, folderWiki)
    cat('New folder created\n')
  }
  
  # Push the file to the current folder
  file <- File(paste0(Path, fileName), parentId = parentId)
  file <- synStore(file)
  cat(fileName, 'pushed\n')
  
  # Add an activity
  rawData <- getEntity(synId)
  used(file)<-list(list(entity = myCode, wasExecuted = TRUE),
                   list(entity = rawData, wasExecuted = FALSE))
  file <- synStore(file)
  
  # Update the folder markdown with the new image to display.
  folderWiki <- synRestGET(folderWikiUri)
  folderUpdateWikiUri <- sprintf("/entity/%s/wiki/%s", parentId, folderWiki$id)
  folderWiki$attachmentFileHandleIds[length(folderWiki$attachmentFileHandleIds)+1] <- file@fileHandle$id
  newMarkdown <-  paste0(paste0(.mkdText(fileName), '\nImage Location: ', propertyValue(file, 'id')),
                         '\n${image?synapseId=', propertyValue(file, 'id'),
                         '&align=None&scale=80}')
  folderWiki$markdown <- paste0(folderWiki$markdown, newMarkdown, '\n')
  
  # now 'push' the wiki to Synapse
  folderWiki <- synRestPUT(folderUpdateWikiUri, folderWiki)
  cat('wiki pushed\n')
  cat('\n')
}

# Helper function: add legends.
.mkdText <- function(fileName){
  if(grepl('Original', fileName)) markdownText <- "###Original PCA:\n*PCA on samples using the full matrix.*"
  if(grepl('PCAfilt', fileName)) markdownText <- "###PCA on features:\n*Using 2nd and 3rd axes, high zero frequencies are located at the center.*"
  if(grepl('samplePCAfilt', fileName)) markdownText <- "###PCA according to filtering:\n*PCA on samples using different levels of filtering (left: using selected probes, right: using rejected probes).*"
  if(grepl('sqDist', fileName)) markdownText <- "###Distances from the center:\n*Distribution of distances (log(d^2)), according to the frequence of zero.*"
  if(grepl('Trace', fileName)) markdownText <- "###Information curve:\n*Samples spread increasing according to the information brought by new probes.*"
  return(markdownText)
}

# Using the cache
# Save in synapse: it works using .rds or .rbin
save(fileName, file = file.path(tempdir(), paste0(fileName, '.rds')))
profileS4 <- File(file.path(tempdir(), paste0(fileName, '.rds')), parentId = 'syn1864121')
profileS4 <- synStore(profileS4)
