
.getValid <- function(x, maxSize = 5){
  tmptab <- table(x)
  commonValue <- as.numeric(names(which(tmptab > maxSize)))
  if(length(commonValue)>0) return(which(x != commonValue))
  else return(seq(1, length(x)))
}

.netScores <- function(Beta){
  Ranks <- lapply(1:ncol(Beta), function(x){nrow(Beta)-rank(abs(Beta[,x]))+1})
  Ranks <- do.call(cbind, Ranks)
  Ranks <- apply(as.data.frame(Ranks), 1, mean, na.rm = T)
  Freq <- apply(Beta, 1, function(x){sum(x!=0)})
  aveBeta <- rowSums(Beta)/ncol(Beta)
  Score <- abs(aveBeta)*Freq
  #Score <- ifelse(Score>1, Score, 1)
  return(list(Beta = Beta, aveBeta = aveBeta, Ranks = Ranks, Freq = Freq, Score = Score))
}
.netModel0 <- function(GE, Resp, B = 100, Nfold = 10, Nlambda = 100){
  output <- lapply(1:B, function(b){
    trainIndex <- sample(1:length(Resp), length(Resp)*0.70)
    trainSet <- GE[ ,trainIndex]
    yTrain <- Resp[trainIndex]
    testSet <- GE[ ,-trainIndex]
    yTest <- Resp[-trainIndex]
    cvnetModel <- cv.glmnet(x = t(trainSet[,bIndex]), y = yTrain[bIndex],
                            nfolds = Nfold, nlambda = Nlambda,
                            family = 'gaussian')
    bestL <- cvnetModel$lambda.min
    bestB <- cvnetModel$glmnet.fit$beta[, which(cvnetModel$lambda == bestL)]
    fit <- predict(cvnetModel, newx = t(testSet[geneIndex,]), s = bestL)
    lmTest <- lm(fit ~ yTest)
    r <- summary(lmTest)$adj.r.squared
    pearsCor <- as.numeric(cor(fit, yTest))
    return(pearsCor)
  })
  output <- do.call(c, output)
  return(output)
}
  
.netModel1 <- function(GE, resp, B = 100, Nfold = 10, Nlambda = 100, fileName,...){
  cat('build model\n')
  #  GE <- ccleGE
  #  resp <- respIC
  #  B = 10
  trainIndex <- sample(1:length(resp), length(resp)*0.70)
  trainSet <- GE[ ,trainIndex]
  yTrain <- resp[trainIndex]
  testSet <- GE[ ,-trainIndex]
  yTest <- resp[-trainIndex]
  
  train <- lapply(seq(1, B), function(x){
    bIndex <- sample(1:ncol(trainSet), replace = TRUE)
    cvnetModel <- cv.glmnet(x = t(trainSet[,bIndex]),
                            y = yTrain[bIndex],
                            nfolds = Nfold, nlambda = Nlambda,
                            family = 'gaussian')
    bestL <- cvnetModel$lambda.min
    bestB <- cvnetModel$glmnet.fit$beta[, which(cvnetModel$lambda == bestL)]
    # Keep the beta values and rank them.
    cat(x, '\tlambda:', bestL, '\tnonZero beta:', sum(bestB != 0), '\tcvm:', min(cvnetModel$cvm),'\n')
    bestB
  })
  
  # Scores features
  trainScores <- .netScores(do.call(cbind, train))
  score <- trainScores$Score
  geneIndex <- which(score > 0)
  cvnetModel <- cv.glmnet(x = t(trainSet[geneIndex,]), y = yTrain, nfolds = 10, nlambda = 100, family = 'gaussian') 
  bestL <- cvnetModel$lambda.min
  bestB <- cvnetModel$glmnet.fit$beta[, which(cvnetModel$lambda == bestL)]
  
  fit <- predict(cvnetModel, newx = t(testSet[geneIndex,]), s = bestL)
  lmTest <- lm(fit ~ yTest)
  r <- summary(lmTest)$adj.r.squared
  pearsCor <- as.numeric(cor(fit, yTest))
  #   png(file = paste0(getwd(), '/',fileName,'.png'), width = 800, height = 800)
  #   par(cex.main= 2, cex.axis = 1.75, cex.lab = 1.5)
  #   plot(sort(yTest), fit[order(yTest)], pch = 19, 
  #        xlab = 'Observed', ylab = 'Fitted', col = 'grey25', asp = 1,...)
  #   abline(0, 1, lty = 3)
  #   abline(lmTest, lwd = 3, col = 'red')
  #   legend('topleft',
  #          legend = c(paste('#Genes:', sum(bestB != 0)),
  #                     paste('adj.r.sq =', round(r, 3)),
  #                     paste('Pearson =', round(pearsCor, 3))),
  #          cex = 1.2, bty = 'n')
  #   dev.off()
  #  return(list(trainScore = trainScores, model = cvnetModel, bestL = bestL, bestB = bestB))
  return(pearsCor)
}

.netModel2 <- function(GE, resp, B = 100, Nfold = 10, Nlambda = 100, kMax = 500, fileName,...){
  cat('build model\n')
  #  GE <- ccleGE
  #  resp <- respIC
  #  B = 10
  trainIndex <- sample(1:length(resp), length(resp)*0.70)
  trainSet <- GE[ ,trainIndex]
  yTrain <- resp[trainIndex]
  testSet <- GE[ ,-trainIndex]
  yTest <- resp[-trainIndex]
  
  train <- lapply(seq(1, B), function(x){
    bIndex <- sample(1:ncol(trainSet), replace = TRUE)
    cvnetModel <- cv.glmnet(x = t(trainSet[,bIndex]),
                            y = yTrain[bIndex],
                            nfolds = Nfold, nlambda = Nlambda,
                            family = 'gaussian')
    bestL <- cvnetModel$lambda.min
    bestB <- cvnetModel$glmnet.fit$beta[, which(cvnetModel$lambda == bestL)]
    # Keep the beta values and rank them.
    cat(x, '\tlambda:', bestL, '\tnonZero beta:', sum(bestB != 0), '\tcvm:', min(cvnetModel$cvm),'\n')
    bestB
  })
  
  # Scores features
  trainScores <- .netScores(do.call(cbind, train))
  
  # Testing Beta matrices
  #Ranks <- trainScores$Ranks
  #ordRanks <- order(Ranks, decreasing = FALSE)
  score <- trainScores$Score
  ordScore <- order(score, decreasing = TRUE)
  #Freq <- trainScores$Freq
  #ordFreq <- order(Freq, decreasing = TRUE)
  
  Rvalues <- c()
  bestK <- bestModel <- NULL
  bestCor <- 0
  for(k in seq(2, min(kMax, sum(ordScore>0)), by = 2)){
    #    if(k%%100 == 0) cat('k:', k, '\n')
    geneIndex <- ordScore[1:k]
    cvnetModel <- cv.glmnet(x = t(trainSet[geneIndex,]), y = yTrain, nfolds = 10, nlambda = 100, family = 'gaussian') 
    bestL <- cvnetModel$lambda.min
    bestB <- cvnetModel$glmnet.fit$beta[, which(cvnetModel$lambda == bestL)]
    if(k%%100 == 0)
      cat('k:', k, ',\tlambda:', bestL, ',\tnum genes:', sum(bestB!=0), '\n')
    #   
    #     #cat('test model\n')
    fit <- predict(cvnetModel, newx = t(trainSet[geneIndex,]))#, s = bestl)
    lmTest <- lm(fit ~ yTrain)
    r <- summary(lmTest)$adj.r.squared
    pearsCor <- as.numeric(cor(fit, yTrain))
    Rvalues <- rbind(Rvalues, cbind(K = k, lambda = bestL, numGene = sum(bestB>0), rsq = r, pearsonCor = pearsCor))
    if(!is.na(pearsCor) & pearsCor > bestCor){
      bestK <- k
      bestCor <- pearsCor
      bestModel <- cvnetModel
    }
  }
  Rvalues <- as.data.frame(Rvalues)
  #    plot(Rvalues$K, Rvalues$V3)
  #    lines(Rvalues$K, Rvalues$V3, col = 'steelblue')
  #    maxC <- Rvalues$K[which.max(Rvalues$pearsonCor)]
  
  #    cvnetModel <- cv.glmnet(x = t(trainSet[ordScore[1:maxC],]), y = yTrain, nfolds = 10, nlambda = 100, family = 'gaussian') 
  bestL <- bestModel$lambda.min
  bestB <- bestModel$glmnet.fit$beta[, which(bestModel$lambda == bestL)]
  fit <- predict(bestModel, newx = t(testSet[ordScore[1:bestK],]))#, s = bestL)
  lmTest <- lm(fit ~ yTest)
  r <- summary(lmTest)$adj.r.squared
  pearsCor <- as.numeric(cor(fit, yTest))
  #     png(file = paste0('/Users/fredcommo/Documents/MyProjects/CellLines/',fileName,'.png'), width = 800, height = 800)
  #     par(cex.main= 2, cex.axis = 1.75, cex.lab = 1.5)
  #     plot(sort(yTest), fit[order(yTest)], pch = 19, 
  #          xlab = 'Observed', ylab = 'Fitted', col = 'grey25', asp = 1,...)
  #     abline(0, 1, lty = 3)
  #     abline(lmTest, lwd = 3, col = 'red')
  #     legend('topleft',
  #            legend = c(paste('#Genes:', sum(bestB != 0)),
  #                       paste('adj.r.sq =', round(r, 3)),
  #                       paste('Pearson =', round(pearsCor, 3))),
  #            cex = 1.2, bty = 'n')
  #     dev.off()
  #  return(list(Rvalues, trainScore = trainScores, model = bestModel, bestL = bestL, bestB = bestB))
  return(pearsCor)
}
