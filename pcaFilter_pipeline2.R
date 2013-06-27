# PCA pipeline2
# compute the frequences of excluded probes
require(VennDiagram)
require(hgu133plus2.db)
require(GO.db)

getIntersect <- function(List, Type){
  idx <- which(tissues == Type)
  probes <- unlist(List)
  for(i in idx)
    probes <- intersect(probes, List[[i]])
  return(probes)
}
#
setwd('/Users/fredcommo/Documents/MyProjects/ProjetACP/PCA_results/')
fileList <- list.files()
fileList <- fileList[grep('.rds$', fileList)]
#
tissues <- lapply(fileList, function(f) unlist(strsplit(f, '_'))[1])
tissues <- do.call(c, tissues)
Names <- lapply(fileList, function(f) unlist(strsplit(f, '_'))[2])
Name <- do.call(c, Names)
#
F02 <- lapply(fileList, function(f){
  cat(f, '\n')
  result <- get(load(f))
  result$rejectList$F0.2$pcaR
})
names(F02) <- Names
#
for(Type in c('Lung', 'Breast', 'CRC')){
  idx = which(tissues == Type)
  n <- length(idx)
  venn.diagram(as.list(F02[idx]), fill = 2:(n+1), alpha = rep(0.1, n),
               cex = 1.5, cat.fontface = 4, lty = 1, fontfamily = 1,
               main = Type, main.cex = 2,
               filename = paste0(Type, '_VennDiagram_F02.png'))
  }
#
lungList <- getIntersect(F02, 'Lung'); length(lungList)
breastList <- getIntersect(F02, 'Breast'); length(breastList)
crcList <- getIntersect(F02, 'CRC'); length(crcList)
venn.diagram(list('Lung\n(n=3)' = lungList, 'Breast\n(n = 4)' = breastList, 'CRC(n = 3)' = crcList),
             fill = 2:4, alpha = rep(0.1, 3),
             cex = 1.5, cat.fontface = 4, lty =1, fontfamily =1,
             filename = paste0('AllTissues_VennDiagram_F02.png'))

crit <- 'var50'
rejectList <- lapply(fileList, function(f){
  cat(f, '\n')
  result <- get(load(f))
  result$rejectList$var50
})
names(rejectList) <- Names
#
for(Type in c('Lung', 'Breast', 'CRC')){
  idx = which(tissues == Type)
  n <- length(idx)
  venn.diagram(as.list(rejectList[idx]), fill = 2:(n+1), alpha = rep(0.1, n),
               cex = 1.5, cat.fontface = 4, lty = 1, fontfamily = 1,
               main = Type, main.cex = 2,
               filename = paste0(Type, '_VennDiagram_', crit,'.png'))
}
#
lungList <- getIntersect(rejectList, 'Lung'); length(lungList)
breastList <- getIntersect(rejectList, 'Breast'); length(breastList)
crcList <- getIntersect(rejectList, 'CRC'); length(crcList)
venn.diagram(list('Lung\n(n=3)' = lungList, 'Breast\n(n = 4)' = breastList, 'CRC(n = 3)' = crcList),
             fill = 2:4, alpha = rep(0.1, 3),
             cex = 1.5, cat.fontface = 4, lty =1, fontfamily =1,
             filename = paste0('AllTissues_VennDiagram_', crit,'.png'))

# Check annotations
lungOnly <- lungList[!is.element(lungList, breastList) & !is.element(lungList, crcList)]
lungOnly <- select(hgu133plus2.db, lungOnly, c("SYMBOL", "GO"), "PROBEID")
lungTerms <- select(GO.db, lungOnly$GO, "TERM", "GOID")
head(lungTerms)

