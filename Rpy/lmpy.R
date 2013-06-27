lmpy <- function(x, y){
	test <- lm(y~x)
	plot(x, y)
	abline(test)
	print (summary(test))
}
