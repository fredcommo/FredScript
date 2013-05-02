#
op <- par(no.readonly = TRUE)
require(survival)
grp <- factor(rep(c('A', 'B'), each = n/2))

n = 100
p = 1
#set.seed(34326)
x <- lapply(1:p, function(x){
  rnorm(n, sample(seq(1, 2, len = 100), 1), sample(seq(0.1, 0.2, len = 10), 1))
  rnorm(n, 3, 0.5)
})
x <- do.call(rbind, x)
Beta1 <- c(seq(-3, -2, len = 100), seq(2, 3, len = 100))
Beta2 <- seq(-2, -3, len = 100)
#LambdaT <- seq(0.02, 0.002, len = 100)
#LambdaC <- seq(0.1, 0.2, len = 100)

x <- matrix(rnorm(n, 2, 0.5),1)
beta1 = 2 #rowMeans(x)/(p) #rep(1.5, p) #sample(Beta1, p) #.5
beta2 = -2*beta1 #sample(Beta2, 1)#-2.5
B <- matrix(c (beta2, beta1), 1)
lambdaT = 0.2 beta1/5 # 0.2 #sample(LambdaT, 1)#.02 								# baseline hazard
lambdaC = 0.2 # beta1/4  								# hazard of censoring
shape = 2
expand = 2.5e-3
#  
B%*%rbind(rep(1, n), x)

#T = rweibull(n, shape=shape, scale=lambdaT*exp(-(beta1*c(0, 1)[grp] + beta2))) + 2.5e-3 	# Event time
#T = rweibull(n, shape=shape, scale=lambdaT*exp(-(beta1*x + beta2))) + expand   # Event time
T = rweibull(n, shape=shape, scale=exp(-(B%*%rbind(rep(1, n), x)))) + expand   # Event time
C = rweibull(n, shape=shape*2, scale=lambdaC)   								#censoring time
Time = pmin(T,C)  													#observed time is min of censored and true
event = (Time==T)*1   												# set to 1 if event is observed
#plot(Time, event)
st <- Surv(Time*1e2, event)
st
km <- survfit(st ~ ifelse(x[1,] < median(x[,1]), 'lo', 'hi'))
plot(km)

#par(mfrow = c(3, 2))
for(i in 1:p){
  km <- survfit(st~ifelse(x[i,]<median(x[i,]), 'lo', 'hi'))
  plot(km, main = mean(x[i,]))
}
par(op)

Data <- as.data.frame(x)
rownames(Data) <- paste0('gene', 1:nrow(x))
colnames(Data) <- paste0('sample', 1:ncol(x))

model = coxph(st ~ ., data = as.data.frame(t(Data)))
summary(model)

par(mfrow = c(4, 2))
for(i in 1:8){
  km <- survfit(st~ifelse(x[i,]<median(x[i,]), 'lo', 'hi'))
  plot(km)
}

require(rbsurv)

#z <- cbind(samples$Age, samples$TumourType)
#z <- cbind(factor(cutClust), samples$TumourType)
fit <- rbsurv(time = Time, status = event, x = Data, 
              gene.ID = rownames(Data), method = "efron",
              max.n.genes = 20, n.iter = 10, nfold = 3)
fit$model

source('/Users/fredcommo/Documents/MyProjects/Fred_Scripts/heatmap.3.R')
Values <- t(x)
heatmap.3(Values,
          scale = "column", Method = "ward",
          col = colorpanel(100, "darkblue", 'grey95', "orange"),
          breaks = seq(-2, 2, len = 101),
          rowsep = 0:(nrow(Values)), colsep = 1:(ncol(Values)-2),
          sepcolor = "grey25", sepwidth=c(0.001,0.001),
          #ColSideColors = probesCol[bestSubGrp],
          #RowSideColors = subGrpCol,
          key = TRUE, keysize = 0.75, trace = "none", density.info = "none")

.generateSurv <- function(Data, shape = 1.2, expand = 2.5e-3,
                          beta1 = 1.35, beta2 = -2.5,
                          lambdaT = 3e-3 lambdaC = 2e-2){
  n <- ncol(Data)
  T = rweibull(n, shape=shape, scale=lambdaT*exp(-(beta1*c(0, 1)[grp] + beta2))) + expand   # Event time
  C = rweibull(n, shape=shape, scale=lambdaC)   								#censoring time
  Time = pmin(T,C)  													#observed time is min of censored and true
  event = (Time==T)*1   												# set to 1 if event is observed
  return
}