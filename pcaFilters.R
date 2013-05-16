####################################################
# Contains the 3 functions to use to perform a PCA-based filtering
# eset is a (p,n) matrix with p genes by rows, and n samples by columns
# pcaGenes <- prcomp(eset)
# score <- pcaTrace1.1(eset, pcaGenes)
# info <- pcaInfo(score)
# select <- pcaSelect(score, p = 0.05) # returns the indices corresponding to 5% of the information
####################################################


pcaTrace1.1 <- function(Data, PCA, Dim = 2:3,...){
  # Data: the original data set
  # PCA: the PCA = prcomp(Data), so a pca on probes.
  # Dim: the dimensions to consider on PCA$x
  # Plot the information curve, and returns a list:
    # m = nrow(Data)
    # n = ncol(Data)
    # PCdim = ncol(X)
    # Dist = D, the vector of distances
    # Score = Score, the Trace scores according to alpha, as a data set
    # lModel = model, the Richard model parameters.
  
  source('/Users/fredcommo/Documents/MyProjects/FredScripts/Richard.w5PL.v2.R')
  X = scale(PCA$x[,Dim])
  X <- as.data.frame(X)
  D <- rowSums(X^2)
  #  a.values <- 10^(seq(log10(1e-4), log10(1e-1), len = 12))
  a = 0.1
  for(i in 1:8) a <- c(a, a[length(a)]*2)
  a.values <- c(1e-4, 1e-3, 1e-2, a)
  
  cat('Testing\n')
  Score <- lapply(a.values, function(a){
    cat('\talpha:', a)
    alpha = 10^(-a)
    Pmax = 1 - alpha
    Q <- qchisq(p = Pmax, df = ncol(X))
    inform <- which(D <= Q)
    lInf <- length(inform)
    
    if(lInf>=2 & lInf<nrow(Data)){
      subData <- Data[inform,]
      acpTest <- prcomp(t(subData))
      tmpTrace <- sum(diag(var(acpTest$x[,1:min(10, ncol(acpTest$x))])))
      tmpScore <- c(aValues = a, nProbes = lInf, Trace = tmpTrace)
      cat('\t#Probes:', lInf, '\tTrace:', tmpTrace, '\n')
      tmpScore
      }
    }
  )
  
  Score <- as.data.frame(do.call(rbind, Score))
  
  # Full matrix
  acpTest <- prcomp(t(Data))
  last.a <- Score$aValues[nrow(Score)]
  tmpTrace <- sum(diag(var(acpTest$x[,1:min(10, ncol(acpTest$x))])))
  tmpScore <- cbind(aValues = c(last.a*2, last.a*4), nProbes = rep(nrow(Data), 2), Trace = rep(tmpTrace, 2))
  Score <- rbind(Score, tmpScore)
  
  cat('\n')
  rownames(Score) <- seq(1, nrow(Score))
  Score <- as.data.frame(Score)
  
  x <- as.numeric(log10(Score$aValues))
  y <- as.numeric(Score$Trace)
  if(any(is.na(y) | is.na(x))){
    na.index <- which(is.na(y) | is.na(x))
    y <- y[-na.index]; x <- x[-na.index]
    }
  y = y/(max(y)*1.01)*100
  if(any(y <=0)) y[y<=0] <- 1e-3
  model <- Richard.w5PL.v2(x, y, w = 0.25, Plot = TRUE, add.points = TRUE,
                           xlab = expression(-Log10(aValues)), ylab = 'Information (%)',...)  
  return(list(m = nrow(Data), n = ncol(Data), PCdim = ncol(X), Dist = D, Score = Score, lModel = model))
}

pcaInfo <- function(pcaScore){
  # pcaScore: pcaTrace output
  # Returns the information table containing the number of informative probes given th proportion of information
  Fmax = pcaScore$lModel$top
  Fb = pcaScore$lModel$bottom
  xmid = pcaScore$lModel$xmid
  b = pcaScore$lModel$scal
  d = pcaScore$lModel$d
  informTable <- c()
  for(i in seq(0.05, 1, by = 0.05)){
    yTarg = (Fmax + Fb)*i
    xTarg = xmid - b*log(((Fmax - Fb)/(yTarg - Fb))^(1/d) - 1)
    
    aTarg <- 10^(xTarg)
    alpha = 10^(-aTarg)
    Pmax = 1 - alpha
    Q <- qchisq(p = Pmax, df = pcaScore$PCdim)  
    inform <- which(pcaScore$Dist >= Q)
    nInform <- length(inform)
    nNonInform <- pcaScore$m - nInform
    informTable <- rbind(informTable,
                         c(Prop = i, Inform = nInform, nonInform = nNonInform, propInform = nInform/pcaScore$m))
  }
  return(as.data.frame(informTable))
}

pcaSelect <- function(pcaScore, p = 0.05){
  # pcaScore: pcaTrace output
  # p: the proprtion of information required
  # Returns the indices of the slected features, according to p.
  Fmax = pcaScore$lModel$top
  Fb = pcaScore$lModel$bottom
  xmid = pcaScore$lModel$xmid
  b = pcaScore$lModel$scal
  d = pcaScore$lModel$d
  
  yTarg = (Fmax + Fb)*p
  xTarg = xmid - b*log(((Fmax - Fb)/(yTarg - Fb))^(1/d) - 1)
  
  aTarg <- 10^(xTarg)
  alpha = 10^(-aTarg)
  Pmax = 1 - alpha
  Q <- qchisq(p = Pmax, df = pcaScore$PCdim)
  inform <- which(pcaScore$Dist >= Q)
  return(inform)
}

