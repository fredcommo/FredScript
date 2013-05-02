# Function
Sophie.test.NewCells <- function(data, Test, Ref, Drug.test=c("SU","ZD"), expandlo = 1.5, expandhi = 3, leg.pos=c(-1.7, -1.9)){
	
	p.prot = as.character(ABlist$Antibody[Test])	
	prot = as.character(ABlist$Antibody[Ref])					
	index.prot <- which(colnames(data)==prot)			# which(fileInfo$Ab==prot)
	index.p.prot <- which(colnames(data)==p.prot)		# which(fileInfo$Ab==prot)

	n <- length(Drug.test)
	par(mfrow=c(1,n))

	# sdata <- split(data[,c(1:5, index.prot, index.p.prot)], data$Replic)
	# R <- length(sdata)
	# All.Res <- c()


	# for (r in 1:R){
	# 	data <- sdata[[r]]

	#	for(D in 1:n){

	index.D <- which(data$Drug== Drug.test | data$Drug== "Placebo")		# Valeurs pour ZD & Placebo: paste("Placebo",Drug.test, sep="_")
	sub_D <- data[index.D, c(1:5, index.prot, index.p.prot)]							#diff2[index.ZD,]
	rownames(sub_D) <- seq(1, nrow(sub_D))
		
	xProt <- log2(sub_D[,which(colnames(sub_D)==prot)])
	xpProt <- log2(sub_D[,which(colnames(sub_D)==p.prot)])
	delta <- xpProt - xProt
	ratio <- cbind.data.frame(sub_D, delta)

	drug <- as.factor(as.character(ratio$Dose))	# ratio$Drug
	dil <- ratio$Dilution

# Part1 : don't consider replicates
	rlm.results <- c()

	for (t in 1:nlevels(data$Time)){
		time=levels(data$Time)[t]
		D_time <- ratio[which(ratio$Time==time), ]
	

		# Estimateur de la valeur moyenne par rlm

		for (d in 1:nlevels(drug)){
			test <- D_time$delta[which(D_time$Dose==levels(drug)[d])]
			dil.test <- log2(D_time$Dilution[which(D_time$Dose==levels(drug)[d])])

			rob <- lm.rob(dil.test, test) # ; rlmB <- rlm(testB~dilB)
			wi <- rob$weights
			s2 <- sum(wi*(predict(rob)-test)^2)/sum(wi)
			s <- sqrt(s2*(1+1/sum(wi)))

			newx = (min(dil.test) + max(dil.test))/2 #; newxB <- (min(dilB) + max(dilB))/2
			estim <-predict(rob, newdata=list(x=newx))

			rlm.results<-rbind(rlm.results, c(time, levels(drug)[d], estim, s))
			}	# end for levels(drug)drug
		}	# end for levels(time)


		rlm.results <- as.data.frame(rlm.results)
		colnames(rlm.results) <- c("Time", "Drug", "estimate", "error")
		rlm.results$estimate <- as.numeric(as.vector(rlm.results$estimate))
		rlm.results$error <- as.numeric(as.vector(rlm.results$error))
		rlm.results <- cbind.data.frame(Comp=paste(p.prot," vs. ",prot, sep=""), rlm.results)
# End part1

# Part2: consider replicates for anova
	replic.results <- c()

	for (t in 1:nlevels(data$Time)){
		time=levels(data$Time)[t]

		for(r in 1:nlevels(data$Replic)){
			rep = levels(data$Replic)[r]
			D_time <- ratio[which(ratio$Time==time & ratio$Replic==rep), ]
	

			# Estimateur de la valeur moyenne par rlm

			for (d in 1:nlevels(drug)){
				test <- D_time$delta[which(D_time$Dose==levels(drug)[d])]
				dil.test <- log2(D_time$Dilution[which(D_time$Dose==levels(drug)[d])])

				rob <- lm.rob(dil.test, test) # ; rlmB <- rlm(testB~dilB)
				wi <- rob$weights
				s2 <- sum(wi*(predict(rob)-test)^2)/sum(wi)
				s <- sqrt(s2*(1+1/sum(wi)))

				newx = (min(dil.test) + max(dil.test))/2 #; newxB <- (min(dilB) + max(dilB))/2
				estim <-predict(rob, newdata=list(x=newx))

				replic.results<-rbind(replic.results, c(time, levels(drug)[d], rep, estim, s))
				}	# end for levels(drug)drug
			}	# End for levels(replicates)
		}	# end for levels(time)


		replic.results <- as.data.frame(replic.results)
		colnames(replic.results) <- c("Time", "Drug", "Replic", "estimate", "error")
		replic.results$estimate <- as.numeric(as.vector(replic.results$estimate))
		replic.results$error <- as.numeric(as.vector(replic.results$error))
		replic.results <- cbind.data.frame(Comp=paste(p.prot," vs. ",prot, sep=""), replic.results)
# End Part2

		# Analyse par Anova

		aov1 <- aov(estimate~Drug*Time, data=replic.results)
		p.drug <- unlist(summary(aov1))[17]
		p.time <- unlist(summary(aov1))[18]

		# Graphe
		rlm.results$Time <- relevel(rlm.results$Time, ref = "4h")
		rlm.results <- rlm.results[order(rlm.results$Time),]
		s.val <- rlm.results$error
		estim <- rlm.results$estimate
		min.y <- min(estim - expandlo*s.val)
		max.y <- max(estim + expandhi*s.val)

		col.palette <- c(rep("indianred1",4), rep("slateblue3",4), rep("violetred",4))	#c(rep("plum",4), rep("mediumorchid1",4), rep("darkviolet",4))
		boxplot(rlm.results$estimate~rlm.results$Drug:rlm.results$Time,  border=NA,
			ylim=range(min.y, max.y), names=rlm.results$Drug, 
			main=unique(rlm.results$Comp),
			cex.main = 2, cex.sub = 2, cex.axis = 1.5, par(bg = "ivory"))

		for (i in 1:nrow(rlm.results)){
			points(i, estim[i], pch=18, cex=4, col=col.palette[i])
			segments(x0=i, x1=i, y0=estim[i]-s.val[i], y1=estim[i]+s.val[i], lwd=3, col=col.palette[i])
			segments(x0=i-0.1, x1=i+0.1, y0=estim[i]-s.val[i], y1=estim[i]-s.val[i], lwd=3, col=col.palette[i])
			segments(x0=i-0.1, x1=i+0.1, y0=estim[i]+s.val[i], y1=estim[i]+s.val[i], lwd=3, col=col.palette[i])
			}

		for (i in seq(1,nrow(rlm.results), by=4))
			lines(estim[i:(i+3)]~seq(i, (i+3)), lwd=3, col=col.palette[i])

		legend1 <- paste("Drug effect: p-value =", signif(p.drug,3))	
		legend2 <- paste("Time effect: p-value =", signif(p.time,3))	

		text(x = c(0.5, 0.5), y = leg.pos, offset = 0, pos = 4, 
				labels=c(legend1, legend2), font=4, cex=1.5, col="royalblue4")
				

		return(list(res =replic.results, Summary=summary(aov1)))

		#All.Res <- rbind(All.Res, rlm.results)
		#}	# end for Drug.test
	#}	# end for replicates test

	#if(length(Drug.test)==1) aov2<-aov(estimate~Dose*Time, data=All.Res)
	#else aov2<-aov(estimate~Drug+Dose+Time, data=All.Res)

	#return(list(All.Res = All.Res, summary(aov2)))
}
