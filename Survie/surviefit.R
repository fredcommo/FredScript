surviefit<-function(time, censor, class, set, thedata, ccol, conftype, Xlim, Xlab, Wlab)
{ts_os <- survfit( Surv(time, censor)~ 1, 
conf.type=conftype,
subset=set, data=thedata)
plot(ts_os, col=ccol, xlim=Xlim,		#c(0, 150), 
xlab=Xlab, ylab=Wlab)

#"Survival Distribution Function")

print(ts_os)}



