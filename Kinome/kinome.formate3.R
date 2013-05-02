kinome.formate3<-function(X, Col = "")

# sort values and build a data.frame in a same format as on the slide

# input : col datas 954 values

# output: 384 rows with well, kinase info, duplicates, mean, sd


{
	Empty.index <- which(X$SpotCoord=="Empty")

	newX <- X[-Empty.index,]
	for (i in 1:6) newX[,i] <- as.character(newX[,i])

	output <- c()

	for (i in seq(1, nrow(newX), by=2))
		output <- rbind(output, c(Rep1 = newX[i, which(colnames(newX)==Col)], Rep2=newX[i+1, which(colnames(newX)==Col)]))

	output <- cbind.data.frame(newX[seq(1, nrow(newX), by=2), 1:8], output)
	output
}

