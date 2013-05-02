
myNet.v3 <- function(Data, Resp, OutFun = "purelin", start = 1, n.max = 8, n.Show = 25, n.rep = 10){

	n.max <- min(n.max, floor((2*nrow(Data)-1)/(ncol(Data) + 2)))
	cat("A maximum of", n.max, "hidden neurons will be evaluated.\n")

	net.err <- c()
	result <- NA
	train.err <- Inf
	for(i in start:n.max){ 	# ceiling(i/2),

	net <- newff(n.neurons = c(ncol(Data), i, ceiling(i/2), 1), learning.rate.global = 1e-2, momentum.global = 0.5,
			error.criterium = "LMS", Stao = NA, hidden.layer = "tansig", output.layer = OutFun, method = "ADAPTgdwm")

		
		tmp.err <- c()
		for(j in 1:n.rep){

			index <- sample(1:nrow(Data), ceiling(nrow(Data)*2/3))
			train.dat <- Data[index,]
			train.resp <- Resp[index]
			test.dat <- Data[-index,]
			test.resp <- Resp[-index]
			
			tmp <- train(net, train.dat, train.resp, error.criterium="LMS", report = F, show.step = 100, n.shows = n.Show)
			fit <- sim(tmp$net, test.dat)
			error <- 1/length(test.resp)*sum((fit - test.resp)^2)
			bic <- log(error) #+ log(i)/2
			# bic <- (-2)*log(error) #+ log(i)/2
			tmp.err <- c(tmp.err, bic)
			cat("Neurons:", i, ", rep:", j, "\n")
			}

		net.err <- c(net.err, mean(tmp.err))
		if(mean(tmp.err) < train.err) {train.err <- mean(tmp.err); result <- tmp}
		}

	return(list(model = result, Err = net.err))
}
