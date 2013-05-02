density.limits<-function(model)

# model : mclustModel() result

{
	n <- dim(model$z)[1]
	n.lim <- length(model$parameters$mean)-1
	means <- model$parameters$mean
	sdev <- model$parameters$variance$sigmasq
	prop <- model$parameters$pro

	frontiers <- rep(0, n.lim)
	if(length(sdev)<(n.lim+1)) sdev <- rep(sdev, n.lim+1)

	for (i in 1:n.lim)
	{
	s2<-1/(n-2)*(n*prop[i]*sdev[i]+n*prop[i+1]*sdev[i+1])
	frontiers[i]<-log(prop[i+1]/prop[i])*s2/(means[i]-means[i+1])+(means[i]+means[i+1])/2
	}
	frontiers
}


density.limits2 <- function(x, param)

# param : parameters from gaussian mixture

{
	n <- length(x)
	n.lim <- length(param)/3
	prop <- param[1:2]
	means <- param[3:4]
	sdev <- param[5:6]

	frontiers <- rep(0, n.lim)

	for (i in 1:n.lim)
	{
	 s2 <- 1/(n-2)*(n*prop[i]*sdev[i]+n*prop[i+1]*sdev[i+1])
	frontiers[i] <- log(prop[i+1]/prop[i])*s2/(means[i]-means[i+1])+(means[i]+means[i+1])/2
	}
	frontiers
}

# used in plot.cloud
