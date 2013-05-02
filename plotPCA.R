# plot PCA
plotPCA <- function(PCA, Axes = 1:2, Cex = 3, l = 2, ...){
# Axis: what PCs to visualize
# Cex: what PC to use for th point sizes
  PCcex <- PCA$x[,Cex]
  mp = min(PCcex)
  Mp = max(PCcex)
  pcex = l*(PCcex-mp)/(Mp - mp) + 0.5
  plot(PCA$x[,Axes], cex = pcex,...)
}