test.fct <- function(data, Test, Ref, Drug.test=c("SU","ZD"), AB=ABlist){
	
	data <- data
	rownames(data) <- seq(1, nrow(data))

	p.prot = as.character(AB$Antibody[Test])	
	prot = as.character(AB$Antibody[Ref])					
	index.prot <- which(colnames(data)==prot)			# which(fileInfo$Ab==prot)
	index.p.prot <- which(colnames(data)==p.prot)		# which(fileInfo$Ab==prot)

	n <- length(Drug.test)
	par(mfrow=c(1,n))

	index.D <- which(data$Drug== Drug.test | data$Drug== "Placebo")		# Valeurs pour ZD & Placebo: paste("Placebo",Drug.test, sep="_")
	sub_D <- data[index.D, c(1:5, index.prot, index.p.prot)]							#diff2[index.ZD,]
		
	# xProt <- log2(sub_D[,which(colnames(sub_D)==prot)])
	# xpProt <- log2(sub_D[,which(colnames(sub_D)==p.prot)])
	# delta <- xpProt - xProt
	# ratio <- cbind.data.frame(sub_D, delta)

	# drug <- as.factor(as.character(ratio$Dose))	# ratio$Drug
	# dil <- ratio$Dilution
	# ratio
	sub_D
}