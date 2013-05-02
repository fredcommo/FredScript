res.cox.sort.adj.feature<-function(
data, time, censor, marqnames, TRT, 
grep_term,
path, res.file, resout.file
						)

# grep_term: term in the logistic regression that should have to be adjusted on its p-value
		# ie. TRT, Biom, TRT:Biom

{
# genenames	: vector of gene symbols corresponding to SNIPS
procs<-c("Bonferroni","Holm","Hochberg","SidakSS","SidakSD","BH","BY")

list_marqueur<-marqnames

N<-ncol(data)						# number of SNIPS

pval<-mute<-naN<-rep(NA, N)

SNPb<-data

#######
# modifiy codings for genotyping data
	
for (s in 1: ncol(data))
{
ind<-which(!is.na(data[,s]) & data[,s] !=  9999999)
leBiom<-data[,s]


# is there an interaction treatment by biomarqueur?

outcox<-coxph(Surv(time[ind], censor[ind]) ~ TRT[ind] + leBiom[ind] + TRT[ind] * leBiom[ind], 
		  method = "efron", data=data)

out<-summary(outcox)

##
# OUTPUT the results that should go into an output file
##

output<-cbind(c("TRT", "Biom", paste("TRT", "Biom", sep=":")),
		  out$conf.int[,-2],
		  out$coef[,5],
      	  rep(list_marqueur[s], 3)
		 )

rownames(output)<-c("TRT", list_marqueur[s], paste("TRT", list_marqueur[s], sep=":"))

# dump the output all togheter in one file
write.table(	output, 
			file		=	paste(path, res.file, sep=''), 
			sep		=	"\t", 
			col.names	=	FALSE, append=TRUE)

}		# end s loop


read.table(paste(path, res.file, sep='')	  , header=F)->resultats

colnames(resultats)<-c("Var", "Term", "HR", "LL", "UL", "p.value", "GeneName")

capture.inter<-which(resultats$Term %in% grep_term)	# grep specific regression term
res<-resultats[capture.inter,]


# correct p-values for FWEAR or FDR using procs
adjpp<-mt.rawp2adjp(as.numeric(as.vector(res[,6])),	procs)
res.adj<-cbind(res[adjpp$index,], adjpp$adjp)

write.table(	res.adj, 
			file		=	paste(path, resout.file, sep=''), 
			sep		=	"\t", 
			row.names	=	FALSE,
			col.names	=	TRUE)
}

# ie. path= "C:/Documents and Settings/a_goubar/Analysis/REMAGUS/Resultats/356SNPs/"

