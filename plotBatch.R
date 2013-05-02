plotBatch <- function(SVD, batch, Col = c(1:nlevels(batch))[batch], Pch = c(1:nlevels(batch))[batch],...){
  require(foreach)
  require(iterators)
  X <- SVD$u[ ,1:2]
  means <- apply(X, 2, function(z){by(z, batch, mean, na.rm = TRUE) })
  PC3 <- SVD$u[,3]
  mp = min(PC3)
  Mp = max(PC3)
  pcex = 2*(PC3-mp)/(Mp - mp) + 0.25
  plot(X, pch = 8, col = Col, cex = pcex)#,...)
  for(b in levels(batch)){
    index = which(batch == b)
    tmp <- X[index,]
    M = means[rownames(means) == b]
    if(length(index) == 1){
      x1 = tmp[1]
      y1 = tmp[2]
      }
    else{
      x1 = tmp[,1]
      y1 = tmp[,2]
      }
    segments(x0 = rep(M[1], length(index)), y0 = rep(M[2], length(index)),
             x1 = x1, y1 = y1, col = Col[index])
    }
  points(means, pch = 18, cex = 1.75, col = unique(Col))
  points(means, pch = 5, cex = 1.75)
  legend('topleft', legend = paste0(levels(batch), ' (n=', table(batch), ')'),
         pch = 8, col = unique(Col), bty = 'n')
}
