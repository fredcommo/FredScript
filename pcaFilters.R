####################################################
# Contains the 3 functions to use to perform a PCA-based filtering
# eset is a (p,n) matrix with p genes by rows, and n samples by columns
# pcaGenes <- prcomp(eset)
# score <- pcaTrace1.1(eset, pcaGenes)
# info <- pcaInfo(score)
# select <- pcaSelect(score, p = 0.05) # returns the indices corresponding to 5% of the information
####################################################

generateRandom <- function(Data, p){
  # Returns a (n,p) matrix of premutation values from p randow Data rows
  output <- lapply(seq(1, p),
                   function(i){ if(i%%ceiling(p/10) == 0) cat(i, '\t')
                                x <- as.numeric(Data[sample(1:nrow(Data),1),])
                                return(sample(x))
                   })
  output <- do.call(rbind, output)
  rownames(output) <- paste0('random', seq(1, nrow(output)))
  colnames(output) <- paste0('sample', seq(1, ncol(output)))
  return(output)
}

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
#      tmpTrace <- sum(diag(var(acpTest$x[,1:min(10, ncol(acpTest$x))])))
      tmpTrace <- sum(acpTest$sdev^2)
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
#  tmpTrace <- sum(diag(var(acpTest$x[,1:min(10, ncol(acpTest$x))])))
  tmpTrace <- sum(acpTest$sdev^2)
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
  model <- .Richard.w5PL.v2(x, y, w = 0.25, Plot = TRUE, add.points = TRUE,
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
    yTarg = Fb + (Fmax - Fb)*i
    xTarg = xmid - b*log(((Fmax*1.01 - Fb)/(yTarg - Fb))^(1/d) - 1)
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

.Richard.w5PL.v2 <- function(x, y, w = 0.25, Plot = F, add.points = F, add.intercept = T, pcol = "royalblue1", add.line = F, lcol = "navy", tan.col = "purple",...){ # Xlim = range(x), Ylim = range(y),
  
  # x  			: x-axis values
  # y				: y-axis values
  # w = 0.25			: weights coefficient
  # Plot = F			: if results have to be visualized
  # add.points = F		: add points on the plot
  # add.intercept = T	: add the intercept point on plot
  # pcol = "royalblue1"	: define points color	
  # add.line = F		: add the tangente line
  # lcol = "navy"		: define the regression color
  # tan.col = "purple"	: define the tangente line color
  # Xlim = range(x)		: the x-axis range
  # Ylim = range(y)		: the y-axis range
  # Title = ""		: to add a plot title
  
  
  x <- as.numeric(x)
  y <- as.numeric(y)
  
  if(any(is.na(y) | is.na(x))){
    na.index <- which(is.na(y) | is.na(x))
    y <- y[-na.index]
    x <- x[-na.index]
  }
  
  # Fonction logistique 5PL
Richard <- function(x, Fb, Fmax, b, c, d){
    y <- Fb + Fmax/(1 + exp(-(x-c)/b))^d
    return(y)
  }
  
  # Fonction sce (somme carr? r?sidus) avec pond?rations
  sce.5P <- function(param, xobs, yobs, w) {
    Fb <- 0 #param[1]
    Fmax <- param[2]
    b <- param[3]
    c <- param[4]
    d <- param[5]
    ytheo <- Richard(xobs, Fb, Fmax, b, c, d)
    sq.res <- (yobs - ytheo)^2
    weights <- 1/sq.res^w
    return(sum(weights*sq.res))
  }
  
  # Fonction sce (somme carr? r?sidus) avec pond?rations
  sce.5P.diag <- function(yobs, ytheo, w) {
    sq.res <- (yobs - ytheo)^2
    weights <- 1/sq.res^w
    return(weights)
  }
  
  # initialisation des parametres
  Fb.ini = min(y)
  Fmax.ini = max(y)	#*1.05
  c.ini = (max(x) + min(x))/2
  z <- (y)/(Fmax.ini - y)
  if (any(abs(z)==Inf)) z[abs(z)==Inf] <- NA
  b.ini = coef(lm(x~log(z)))[2]
  # b.ini = 1					
  d.ini = 1
  init <- c(Fb.ini, Fmax.ini, b.ini, c.ini, d.ini)
  
  # Estimation du modele
  best<-nlm(f = sce.5P, p = init, xobs = x, yobs = y, w = w)
  
  # R?cup?ration des param?tres
  best.Fb <- best$estimate[1]
  best.Fmax <- best$estimate[2]
  best.b<-best$estimate[3]
  best.c <- best$estimate[4]
  best.d <- best$estimate[5]
  
  # Diagnostic de r?gression
  yfit <- Richard(x, best.Fb, best.Fmax, best.b, best.c, best.d)
  weights <- sce.5P.diag(y, yfit, w)
  lm.test <- lm(yfit~y)	#, weights = weights)
  r.sq <- summary(lm.test)$adj.r.squared
  lm.slope <-coef(lm.test)[2]
  p.slope <- summary(lm.test)$coefficients[2,4]
  
  # Estimation des valeurs pour graphique
  newx <- seq(min(x), max(x), length=100)						
  newy <- Richard(newx, best.Fb, best.Fmax, best.b, best.c, best.d)
  
  # coordonn?es du pt d'inflexion
  Xflex = best.c + best.b*log(best.d)
  Yflex = best.Fb + best.Fmax*(best.d/(1 + best.d))^(best.d)
  
  # coordonn?es du pt rep = 0.5
  Y50 = (max(newy) + min(newy))/2
  X50 = best.c - best.b*log((best.Fmax/(Y50 - best.Fb))^(1/best.d) - 1)
  
  # pente au pt d'inflexion
  B = best.Fmax/best.b*(best.d/(1 + best.d))^(best.d + 1); B
  A = Yflex  - B*(Xflex); A
  
  # pente finale
  x.ini <- x[1]
  x.end <- x[length(x)]
  y.ini <- Richard(x.ini, best.Fb, best.Fmax, best.b, best.c, best.d)
  y.end <- Richard(x.end, best.Fb, best.Fmax, best.b, best.c, best.d)
  Bf = (best.d*best.Fmax/best.b)*exp(-1/best.b*(x.end-best.c))*(1+exp(-1/best.b*(x.end-best.c)))^(-best.d-1)
  Af = y.end - Bf*x.end
  
  # x.intercept
  Cy0 = -A/B
  
  # Repr?sentations graphiques
  if(Plot){
    plot(y~x, type = 'n',...)
    lines(newy~newx, col = pcol,...)
    lines(I(A+B*x)~x, col = tan.col, lwd = 2)
    abline(h = 0, lwd = 1, lty = 3, col = 'grey75')
    if(add.intercept) points(-A/B, 0, pch = 19, col = "red")
  }
  
  if(add.points) points(y~x, col = pcol,...)
  
  if(add.line){
    lines(newy~newx, col = lcol, lwd = 2)
    lines(I(A+B*x)~x, col = tan.col, lwd = 1)
    if(add.intercept) points(-A/B, 0, pch = 19, col = "red")
  }
  
  return(list(bottom = best.Fb, top = best.Fmax, xmid = best.c, scal = best.b, d = best.d,
              Xflex = Xflex, Yflex = Yflex, slope = B, x.intercept = Cy0, Yini = y.ini, Yend = y.end, end.slope = Bf,
              lm.rsq = r.sq, lm.slope = lm.slope, p.value = p.slope, xfit = newx, yfit = newy))
}

