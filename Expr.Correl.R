Expr.Correl<-function(x,data,method="spearman")

# x 		: reference (expression to be compare to)
# data	: data frame of all the values to compare to the reference


{
	output<-data.frame(Probe=rownames(data),estimate=0,p.value=0)
		for (i in 1:nrow(data)){
			x<-as.numeric(x)+rnorm(length(x),mean=0,sd=1e-10)
			y<-as.numeric(data[i,])+rnorm(length(x),mean=0,sd=1e-10)

			test<-cor.test(x,y,method="spearman",exact=TRUE)
			output$estimate[i]<-test$estimate
			output$p.value[i]<-test$p.value
		}
	output<-output[order(abs(output$estimate),decreasing=TRUE),]

	Pos.Cor<-output[which(output$estimate>0),]
	Pos.Cor<-Pos.Cor[order(Pos.Cor$estimate,decreasing=TRUE),]

	Neg.Cor<-output[which(output$estimate<0),]
	Neg.Cor<-Neg.Cor[order(Neg.Cor$estimate),]


	list(Pos.Cor=Pos.Cor, Neg.Cor=Neg.Cor)
}

library(annaffy)

probeId<-Neg.Cor$Probe[1:10]
probeId<-as.character(probeId)
values<-Neg.Cor$estimate[1:10]

ann.col<-aaf.handler()[-c(4:6)]

ann.table<-aafTableAnn(probeId, "hgu133a", ann.col)
saveHTML(ann.table, "Correlations.html", title="Best correlations (Spearman)")





