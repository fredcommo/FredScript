#
set.seed(34326)
n = 100
grp <- factor(rep(c('A', 'B'), each = n/2))

beta1 = 1.35; beta2 = -2.5
lambdaT = .003 								# baseline hazard
lambdaC = .02  								# hazard of censoring
s = 1.2

T = rweibull(n, shape=s, scale=lambdaT*exp(-(beta1*c(0, 1)[grp] + beta2))) + 2.5e-3 	# Event time
C = rweibull(n, shape=s, scale=lambdaC)   								#censoring time
Time = pmin(T,C)  													#observed time is min of censored and true
event = (Time==T)*1   												# set to 1 if event is observed
plot(Time, event)

st <- Surv(Time*2e3, event)
model = coxph(st ~ grp)
summary(model)
km <- survfit(st~grp)
plot(km)

Data <- cbind.data.frame(Sample = paste0('sample', seq(1, length(event))), Grp = grp, fUp = Time*3e3, Event = event*1)

