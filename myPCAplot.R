myPCAplot <- function(X, Fact, Col = NULL, Expand = 0.5, Leg = TRUE,...){

# X : matrix of new components containing at least 3 dimensions ( <- prcomp()$x)
# Fact : a vector of classes.
# Expand : a constant value to adjust the size of points
# ... : additional graphics parameters

	if(!is.factor(Fact)) Fact <- as.factor(Fact)	
	PC3 <- X[,3]
	mp = min(PC3)
	Mp = max(PC3)
	pcex = 2*(PC3-mp)/(Mp - mp) + Expand
	Cols <- rep(NA, length(pcex))
	ng <- nlevels(Fact)
	LegCol <- c()
	if(is.null(Col)) RGB <- c()
		for(i in 1:ng){
			grpIndex <- which(Fact == levels(Fact)[i])
			
			if(is.null(Col)){
				rgb1 <- (ng-(ng-i+1))/ng
				rgb2 <- sample(seq(0, 1, by = 0.05), 1)
				rgb3 <- 1- rgb1
				RGB <- rbind(RGB, c(rgb1, rgb2, rgb3))
				}
			else rgb1 <- RGB[i, 1]; rgb2 <- RGB[i, 2]; rgb3 <- RGB[i, 3]
			
			Cols[grpIndex] <- rgb(rgb1, rgb2, rgb3, alpha = pcex[grpIndex]/max(pcex))
			LegCol <- c(LegCol, rgb(rgb1, rgb2, rgb3, alpha = 1))
		}
	plot(X[,1:2], pch = 19,
		cex = pcex,
		col = Cols,
		asp = 1,...)
		if(Leg){
			legend("topright", legend = levels(Fact), pch = 19, col = LegCol, bty = "n")
		}
		return(RGB)
}

############################
# Examples

# data(iris)
# pca <- prcomp(iris[,-5])$x
# grp <- iris[,5]
# myPCAplot(pca, grp, main = "Iris example")

# data(crabs)
# values <- crabs[,4:8]
# values <- values[,-3]/values[,3]
# grp <- crabs[,1]:crabs[,2]
# pca <- prcomp(values)$x
# myPCAplot(pca, grp, Expand = 0.9, main = "Crabs example")

