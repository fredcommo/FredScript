setwd('/Volumes/FREECOM HDD/Borrensen_Data/ReAnalyzed_aCGH_UIIsamples_GSE20394')

ftpTable <- read.table('GSE20394_sample_list.txt', header = TRUE, sep = '\t', stringsAsFactors=FALSE)
ftpList <- ftpTable$Sample_supplementary_file
destName <- sapply(ftpList, function(x){tmp <- unlist(strsplit(x, '\\/'))
                                        tmp[length(tmp)]})

setwd('/Volumes/FREECOM HDD/Borrensen_Data/ReAnalyzed_aCGH_UIIsamples_GSE20394/GSE20394_RawData')
for(i in 1:length(ftpList)){
  fileUrl <- ftpList[i]
  download.file(fileUrl,destfile=as.character(destName[i]),method="curl")
  cat(as.character(destName[i]), 'downloaded\n')
}

tar('./GSE20394_RAW.tar', files = list.files())

# MicMa_WZ data GSE19425
setwd('/Volumes/FREECOM HDD/Borrensen_Data/MicMA_WZ')
dir.create('MicMa_RawData')
ftpTable <- read.table('MicMa_filelist.txt', header = TRUE, sep = '\t', stringsAsFactors=FALSE)
fileNames <- ftpTable$Name[-c(1:3)]
fileList <- ftpTable$ftp[-c(1:3)]

setwd('/Volumes/FREECOM HDD/Borrensen_Data/MicMA_WZ/MicMa_RawData/')
for(i in 1:length(fileList)){
  cat('\n', fileNames[i], '...\n')
  fileUrl <- fileList[i]
  download.file(fileUrl, destfile = fileNames[i], method = 'auto')
}


# NCI-60 Gene Expr Raw
setwd('/Users/fredcommo/Documents/Sage/Projects/NCI_60 Cell Lines/NCI60 GeneExpr')
dir.create('MicMa_RawData')
ftpTable <- read.table('GSE22821_AnnotTable.txt', header = TRUE, sep = '\t', stringsAsFactors=FALSE)
fileNames <- ftpTable$Sample_geo_accession
fileList <- ftpTable$Sample_supplementary_file

setwd('/Volumes/FREECOM HDD/Borrensen_Data/MicMA_WZ/MicMa_RawData/')
for(i in 1:length(fileList)){
  cat('\n', fileNames[i], '...\n')
  fileUrl <- fileList[i]
  download.file(fileUrl, destfile = fileNames[i], method = 'auto')
}
