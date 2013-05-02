brackets <- function(xleft, xright, y, hi = 0.25, len = 0, Lty = 3, Lwd = 1, Col="black"){
	segments(x0 = xleft, x1 = xleft, y0 = y - hi, y1 = y + hi, lty = Lty, lwd = Lwd, col = Col)
	segments(x0 = xleft, x1 = xleft + len, y0 = y - hi, y1 = y - hi, lty = Lty, lwd = Lwd, col = Col)
	segments(x0 = xleft, x1 = xleft + len, y0 = y + hi, y1 = y + hi, lty = Lty, lwd = Lwd, col = Col)

	segments(x0 = xright, x1 = xright, y0 = y - hi, y1 = y + hi, lty = Lty, lwd = Lwd, col = Col)
	segments(x0 = xright - len, x1 = xright, y0 = y - hi, y1 = y - hi, lty = Lty, lwd = Lwd, col = Col)
	segments(x0 = xright - len, x1 = xright, y0 = y + hi, y1 = y + hi, lty = Lty, lwd = Lwd, col = Col)
	}
