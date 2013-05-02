VennCounts3<-function (results, include = "both") 
{
    	x <- as.matrix(results)
    	include <- match.arg(include, c("both", "up", "down"))
    	#include<-"both"
	
	x <- sign(switch(include, both = abs(x), up = x > 0, down = x < 0))
    	nprobes <- nrow(x)
    	ncontrasts <- ncol(x)
    	names <- colnames(x)

    	if (is.null(names)) 
        	names <- paste("Group", 1:ncontrasts)

    	noutcomes <- 2^ncontrasts
    	outcomes <- matrix(0, noutcomes, ncontrasts)
    	colnames(outcomes) <- names

    	for (j in 1:ncontrasts) 
		outcomes[, j] <- rep(0:1, times = 2^(j-1), each = 2^(ncontrasts-j))

    	xlist <- list()

    	for (i in 1:ncontrasts) 
		xlist[[i]] <- factor(x[, ncontrasts-i+1], levels = c(0, 1))

    # Extract genelists from Venn zones
	which(xlist[[1]]==0 & xlist[[2]]==0 & xlist[[3]]==1)->list1
	which(xlist[[1]]==0 & xlist[[2]]==1 & xlist[[3]]==0)->list2
	which(xlist[[1]]==1 & xlist[[2]]==0 & xlist[[3]]==0)->list3

    															#list(grpA=list1,grpB=list2,grpC=list3)->discrim

	counts <- as.vector(table(xlist))
    	list(Table=cbind(outcomes, Counts = counts),grpA=list1,grpB=list2,grpC=list3 )	#Grps=discrim
}

