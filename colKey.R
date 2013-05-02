

# Script: Add a color key
colKey <- function(cStart = "white", cMid = NA, cEnd = "purple", n = 1000, Range = c(0, 1), xlab = "Values", ylab = NA,...){
	require(gplots)
	op <- par(no.readonly = T)

	Col <- colorpanel(n, cStart, cMid, cEnd)
	if(is.na(cMid)) Col <- colorpanel(n, cStart, cEnd)
	x = seq(Range[1], Range[2], len = n)
	y = 1

	par(bty = "n", mar = c(15, 10, 15, 10))
		image(x, y, z = matrix(x, n, 1), col = Col, axes = F, xlab = xlab, ylab = ylab,...)
		myTicks <- seq(Range[1], Range[2], len = 5)
	axis(side = 1, at = myTicks, labels = myTicks)
	par(op)
	}

# ex : colKey("yellow", "white", "blue")
# ex : colKey("green3", "black", "red3", n = 50, Range = c(-5, 5), xlab = "Z-values")

