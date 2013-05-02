kinome.formate2<-function(X)

# sort values and build a data.frame in a same format as on the slide

# input : col datas 954 values

# output: 384 rows with well, kinase info, duplicates, mean, sd


{
	Empty.index <- which(X$SpotCoord=="Empty")

	newX <- X[-Empty.index,]
	for (i in 1:6) newX[,i] <- as.character(newX[,i])

	output <- c()

	for (i in seq(1, nrow(newX), by=2))
		output <- rbind(output, c(Rep1 = newX$F635[i], Rep2=newX$F635[i+1]))

	output <- cbind.data.frame(newX[seq(1, nrow(newX), by=2), 1:8], output)
	output
}

