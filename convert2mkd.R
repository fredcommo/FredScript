
# Convert a data.frame into markdown

convert2mkd <- function(myTab, headers = TRUE){
  output = ''
  if (headers){
    head.sep = ''
    foreach(n = iter(names(myTab))) %do% {
    	output = paste(output, n, sep = '|')
      	head.sep = paste(head.sep, '-', sep = '|')
    	}
    output = paste0(output, '|\n', head.sep, '|\n')		
  }
  for (r in 1:nrow(myTab)){
    tmp <- myTab[r,]
    foreach(s=iter(tmp)) %do% {output = paste(output, s, sep = '|')}
    output = paste0(output, '|\n')
 	}
  return(output)
}
