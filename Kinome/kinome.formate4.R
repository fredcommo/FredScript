kinome.formate4<-function(X, Info){

# sort values and build a data.frame in a same format as on the slide
# input : col datas 954 values
# output: 384 rows with well, kinase info, duplicates, mean, sd

formate <- function(Y){
	output <- c()
	for (i in seq(1, length(Y), by=2)) output <- rbind(output, c(Y[i], Y[i+1]))
	output <- apply(output, 1, mean, na.rm=T)
	return(output)	
	}

	newX <- c()

	for (i in seq(1, ncol(X))) newX <- cbind(newX, formate(X[,i]))

	colnames(newX) <- colnames(X)
	newInfo <- Info[seq(1, nrow(X), by=2),]

	return(list(keset.rep = newX, k.info.rep = newInfo))
}

