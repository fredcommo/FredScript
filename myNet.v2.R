myNet.v2 <- function(Data, Resp, n.max = 6, n.Show = 25, n.rep = 10){

	net.err <- c()
	# plot(rep(0, n.runs)~seq(1, n.runs), ylim = range(0, 0.15),type = "n")
	for(i in 1:n.max){
		err <- c()
		for(j in 1:n.rep){
			index <- sample(1:nrow(Data), ceiling(nrow(Data)*2/3))
			train.dat <- Data[index,]
			train.resp <- Resp[index]
			test.dat <- Data[-index,]
			test.resp <- Resp[-index]
			net <- newff(n.neurons=c(ncol(Data), i, ceiling(i/2), 1), learning.rate.global=1e-2, momentum.global=0.5,
				error.criterium="LMS", Stao=NA, hidden.layer="tansig", output.layer="purelin", method="ADAPTgdwm")
			result <- train(net, train.dat, train.resp, error.criterium="LMS", report = F, show.step = 100, n.shows = n.Show)
			fit <- sim(result$net, test.dat)
			error <- 1/length(test.resp)*sum((fit - test.resp)^2)
			bic <- (-2)*log(error) + log(i)/2
			# bic <- (-2)*log(1/sum(Merrors)) + i*log(ncol(Data))
			err <- c(err, bic)
			cat("Neurons:", i, ", rep:", j, "\n")
			}
		net.err <- c(net.err, mean(err))
		}
# legend("topright", legend = paste("nnc =", seq(1, ncc.max)), lty = 1, col = 1:ncc.max)

ncc <- which.min(net.err)
net <- newff(n.neurons=c(ncol(Data), ncc, ceiling(ncc/2), 1), learning.rate.global=1e-2, momentum.global=0.5,
		error.criterium="LMS", Stao=NA, hidden.layer="tansig", output.layer="purelin", method="ADAPTgdwm")
err <- Inf
result <- NA
for(j in 1:n.rep){
	index <- sample(1:nrow(Data), ceiling(nrow(Data)*2/3))
	train.dat <- Data[index,]
	train.resp <- Resp[index]
	test.dat <- Data[-index,]
	test.resp <- Resp[-index]
	net <- newff(n.neurons=c(ncol(train.dat), i, ceiling(i/2), 1), learning.rate.global=1e-2, momentum.global=0.5,
		error.criterium="LMS", Stao=NA, hidden.layer="tansig", output.layer="purelin", method="ADAPTgdwm")
	tmp <- train(net, train.dat, train.resp, error.criterium="LMS", report = F, show.step = 100, n.shows = n.Show)
	fit <- sim(tmp$net, test.dat)
	error <- 1/length(test.resp)*sum((fit - test.resp)^2)
	if(error < err) {err <- error; result <- tmp}
	cat("rep:", j, ", error:", err, "\n")
	}

	return(result)
}
