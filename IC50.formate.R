IC50.formate <- function(dataset, nrep = 3){

	cell.lines <- dataset$CellLine
	doses <- dataset[,2]
	values <- dataset[,-c(1:2)]
	resp <- (values[,1:nrep]-dataset$T0)/(dataset$Ctrl-dataset$T0)
	mean.resp <- apply(resp, 1, function(x){x[x < 0] <- 0.01; prod(x, na.rm = T)^(1/length(x))})
	sd.resp <- apply(resp, 1, sd, na.rm = T)
	nbval <- apply(resp, 1, length)
	SEM <- sd.resp/sqrt(nbval)

	newdataset <- cbind.data.frame(cell.lines, doses, resp, mean.resp, SEM)
	colnames(newdataset) <- c("CellLine", "Dose", paste("Replic", seq(1, nrep), sep=""), "Mean", "SEM")
	return(newdataset)
}
