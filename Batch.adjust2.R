Batch.adjust2 <- function(Eset, info, qnorm = F, Center = T, Scale = T){

# Eset :	matrix(n,m) of expr. values
# info : 	a 2 columns matrix (or data.frame). col#1 = exp ids, col#2 = batch

require(preprocessCore)

	verif <- ifelse(colnames(Eset)==info[,1], "ok", "error")
	errors <- which(verif=="error")
	if(length(errors)>0) stop("data and patients are in different order. Please fix it and try again!")

	cols <- colnames(Eset)
	rows <- rownames(Eset)

	if(qnorm) Eset <- normalize.quantiles(as.matrix(Eset))

	sEset <- apply(Eset, 2, function(x)((x - mean(x, na.rm = T))/sd(x, na.rm = T)))
	batch <- unique(info[,2]); cat("Found", nlevels(factor(batch)), "batche(s)\n")
	nEset <- matrix(0, nrow(Eset), ncol(Eset))
	for (i in 1:length(batch)){
		tmp <- sEset[,info[,2]==batch[i]]
		tmp <- apply(tmp, 1, function(x)((x - mean(x, na.rm = T))/sd(x, na.rm = T)))
		nEset[,info[,2]==batch[i]] <- t(tmp)
		}
	colnames(nEset) <- cols
	rownames(nEset) <- rows
	cat("Your data set have been standardized ;)\n")
	return(nEset)
}

