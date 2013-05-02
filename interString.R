
interString <- function(x, y){

	s1 <- rep(NA, nchar(x))
	s2 <- rep(NA, nchar(y))

	for (i in 1:nchar(x)) s1[i]<- substr(x, start=i, stop=i)
	for (i in 1:nchar(y)) s2[i]<- substr(y, start=i, stop=i)

	for (i in 1:length(s1)){
		c <- s1[i]
		pos <- which(s2==c)
		
		if(length(pos)>0)
			if(s2[pos+1]==s1[i+1] & !is.na(s1[i+j])){
			common <- c(s1[i], s1[i+1])
			j=2
			while(s2[pos+j]==s1[i+j] & !is.na(s1[i+j])){
				common <- c(common, s1[i+j])
				j=j+1
			}
		}
	}
	pept <- c()
	for (i in 1:length(common)) pept <- paste(pept, common[i], sep="")
	pept
}