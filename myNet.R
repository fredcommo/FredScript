myNet <- function(Data, Resp, ncc.max = 8, n.run = 25, n.rep = 10){

err <- c()
# plot(rep(0, n.runs)~seq(1, n.runs), ylim = range(0, 0.15),type = "n")
for(i in 1:ncc.max){
	cat("Estimating...\n")
	for(j in 1:n.rep){
	net <- newff(n.neurons = c(ncol(Data), i, ceiling(i/2), 1), learning.rate.global=1e-2, momentum.global=0.5,
		error.criterium="LMS", Stao=NA, hidden.layer="tansig", output.layer="purelin", method="ADAPTgdwm")
	result <- train(net, Data, Resp, error.criterium="LMS", report = T, show.step = 100, n.shows = n.runs)
	Merrors <- result$Merror
	# lines(Merrors, lwd = 2, col = i)
	bic <- (-2)*log(1/sum(Merrors)) + i
	# bic <- (-2)*log(1/sum(Merrors)) + i*log(ncol(Data))
	err <- c(err, bic)
	cat("model", i, ": trial", j, "\n")
	}
}
# legend("topright", legend = paste("nnc =", seq(1, ncc.max)), lty = 1, col = 1:ncc.max)

ncc <- ceiling(which.min(err)/n.rep); ncc
cat("Best model estimated with", ncc, "neurons\n")
net <- newff(n.neurons=c(ncol(Data), ncc, ceiling(ncc/2), 1), learning.rate.global=1e-2, momentum.global=0.5,
		error.criterium="LMS", Stao=NA, hidden.layer="tansig", output.layer="purelin", method="ADAPTgdwm")
err <- Inf
result <- NA

cat("Training...\n")
for(j in 1:n.rep){	
	tmp <- train(net, Data, Resp, error.criterium = "LMS", report = T, show.step = 100, n.shows = n.runs)
	tmp.err <- (-2)*log(1/sum(tmp$Merror))
	if(tmp.err < err) {err <- tmp.err; result <- tmp}
	cat("rep:", j, "\n")
	}	

# y <- sim(result$net, values)
cat("Done\n")
# return(list(y = y, Archi = ncc))
return(result)
}
