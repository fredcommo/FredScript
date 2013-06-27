require(synapseClient)

hapmap1  <- synGet('syn1929397')

celPath  <- paste0(gsub('\\?.+', '', hapmap1@filePath), 'CEL')
untar(hapmap1@filePath, exdir = celPath)

require(aroma)
require(pd.genomewidesnp.6)
setwd(celPath)
rawData  <- ReadAffy()

# Normalization
eset.mas5  <- mas5(rawData)
eset  <- log2(exprs(eset.mas5))

