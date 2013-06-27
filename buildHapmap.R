# Build hapmap3

require(aroma.affymetrix)
require(oligo)
require(synapseClient)
require(affy)
require(pd.genomewidesnp.6)

query <- synapseQuery("select name, id from entity where parentId == 'syn1929395'")
synIds <- query$entity.id[grepl('syn',query$entity.name)]

celPath  <- '/Users/fredcommo/.synapseCache/HP'
dir.create(celPath)
dir.create(outputDir, recursive = TRUE)

lapply(synIds, function(id){
	hp  <- synGet(id)
	untar(hp@filePath, exdir = celPath)
})

setwd(celPath)
outputDir <- file.path(getwd(),'crlmmResults')
crlmm(list.celfiles(celPath, full.names=TRUE), outputDir)
crlmmOut <- getCrlmmSummaries(outputDir)

rawData  <- ReadAffy(list.celfiles())
eset.mas5  <- mas5(rawData)

