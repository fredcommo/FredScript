require(glmnet)

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

.netModel1 <- function(GE, resp, B = 100, Nfold = 10, Nlambda = 100, alpha = 0.1, fileName,...){
  cat('build model\n')
  #GE <- Data
  #resp <- mekResp
  #B = 20
  if(any(is.na(resp))){
    NAs <- which(is.na(resp))
    GE <- GE[, -NAs]
    resp = resp[-NAs]
  }
  trainIndex <- sample(1:length(resp), length(resp)*0.70)
  trainSet <- GE[ ,trainIndex]
  yTrain <- resp[trainIndex]
  testSet <- GE[ ,-trainIndex]
  yTest <- resp[-trainIndex]
  
  train <- lapply(seq(1, B), function(b){
    bIndex <- sample(1:ncol(trainSet), replace = TRUE)
    cvnetModel <- cv.glmnet(x = t(trainSet[,bIndex]),
                            y = yTrain[bIndex],
                            alpha = alpha, nfolds = Nfold, nlambda = Nlambda,
                            family = 'gaussian')
    bestL <- cvnetModel$lambda.min
    bestB <- cvnetModel$glmnet.fit$beta[, which(cvnetModel$lambda == bestL)]
    # Keep the beta values and rank them.
    cat(b, '\tlambda:', bestL, '\tnonZero beta:', sum(bestB != 0), '\tcvm:', min(cvnetModel$cvm),'\n')
    bestB
  })
  
  # Scores features
  trainScores <- .netScores(do.call(cbind, train))
  score <- trainScores$Score
  geneIndex <- which(score > 0)
  cat(length(geneIndex), 'found with non Zero beta values\n')
  cvnetModel <- cv.glmnet(x = t(trainSet[geneIndex,]), y = yTrain, alpha = alpha, nfolds = Nfold, nlambda = Nlambda, family = 'gaussian') 
  bestL <- cvnetModel$lambda.min
  bestB <- cvnetModel$glmnet.fit$beta[, which(cvnetModel$lambda == bestL)]

  fit <- predict(cvnetModel, newx = t(testSet[geneIndex,]), s = bestL)
  lmTest <- lm(fit ~ yTest)
  r <- summary(lmTest)$adj.r.squared
  pearsCor <- as.numeric(cor(fit, yTest))
  return(list(trainScore = trainScores, model = cvnetModel, bestL = bestL, bestB = bestB,
             boostGenes = geneIndex, y = yTest, yFit = fit, pearsCor = pearsCor))
}

.netModel2 <- function(GE, resp, B = 100, Nfold = 10, Nlambda = 100, kMax = 500, fileName,...){
  cat('build model\n')
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
#  return(list(Rvalues, trainScore = trainScores, model = bestModel, bestL = bestL, bestB = bestB))
  return(pearsCor)
}


