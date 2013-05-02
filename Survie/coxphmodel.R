coxphmodel<-function(time, censor, class, covariables, thedata)
{
outcox<-coxph(Surv(time, censor) ~ class , method = "efron", data=thedata)	
out<-summary(outcox)
cbind(out$conf.int, out$coef[,5])
}

