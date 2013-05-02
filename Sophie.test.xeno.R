# Function
Sophie.test.xeno <- function(data, Test, Ref, Drug.test=c("SU","ZD")){
	n <- length(Drug.test)
	par(mfrow=c(1,n))
	All.Res <- c()

	for(D in 1:n){
		index.D <- which(Drug== Drug.test[D] | Drug== paste("Placebo",Drug.test[D], sep="_"))	# Valeurs pour ZD & Placebo
		sub_D <- data[index.D,]	#diff2[index.ZD,]

		# testé = forme phospho
		# référence = forme totale

		test = Test
		ref = Ref
		graph.range=NA

		p.prot = as.character(ABlist$Antibody[test])	
		prot = as.character(ABlist$Antibody[ref])					

		index.prot <- which(colnames(sub_D)==prot)		# which(fileInfo$Ab==prot)
		index.p.prot <- which(colnames(sub_D)==p.prot)		# which(fileInfo$Ab==prot)


		rlm.results <- c()
			for (t in c(1,2)){
				time=levels(Time)[t]

				D_time <- sub_D[which(sub_D$Time==time), ]

				xA <- log2(D_time[,index.prot])			# toutes les valeurs 'forme totale' au temps T, toute concentration, tout replicat
				xB <- log2(D_time[,index.p.prot])			# toutes les valeurs 'forme phospho' au temps T, toute concentration, tout replicat

				ratio <- data.frame()

				for (i in seq(1,nrow(D_time), by=3)){
					A <- xA[i:(i+2)]
					B <- xB[i:(i+2)]
					D1 <- diag(rep(A, each=3))
					D2 <- diag(rep(B, 3))
					delta <- diag(D2-D1)

					tmp <- data.frame()
					for (j in 1:9) tmp <- rbind.data.frame(tmp, D_time[i,1:4])	
			
					tmp <- cbind.data.frame(tmp, delta)
					ratio <- rbind.data.frame(ratio, tmp)	# tous les delta phospho vs. total, pour toute concentration, au temps T, inclus ZD & Placebo 
					}

				drug <- as.factor(as.character(ratio$Drug))
				# dose <- as.factor(as.character(ratio$Dose))
				# dd <- as.factor(as.character(drug:dose))
				dil <- ratio$Dilution


				# Estimateur de la valeur moyenne par rlm

				for (d in 1:nlevels(drug)){
				#	d.dose <- c("0,1", "0,5", "1")
				#	if (levels(drug)[d]=="Placebo") d.dose <- "0"
		
				#		for (dos in 1:length(d.dose)){
							test <- ratio$delta[which(drug==levels(drug)[d])]
							dil.test <- log2(dil[which(drug==levels(drug)[d])])

							rob <- lm.rob(dil.test, test) # ; rlmB <- rlm(testB~dilB)

							wi <- rob$weights
							s2 <- sum(wi*(predict(rob)-test)^2)/sum(wi)
							s <- sqrt(s2*(1+1/sum(wi)))

							newx = (min(dil.test) + max(dil.test))/2 #; newxB <- (min(dilB) + max(dilB))/2
							estim <-predict(rob, newdata=(x=newx))


							rlm.results<-rbind(rlm.results, c(time, levels(drug)[d], estim, s))
				#			}
					}
				}

		rlm.results <- as.data.frame(rlm.results)
		colnames(rlm.results) <- c("time", "Drug", "estimate", "error")
		rlm.results$estimate <- as.numeric(as.vector(rlm.results$estimate))
		rlm.results$error <- as.numeric(as.vector(rlm.results$error))
		rlm.results <- cbind.data.frame(Comp=paste(p.prot," vs. ",prot, sep=""), rlm.results)

		# Analyse par Anova
		#estimate <- as.numeric(as.vector(rlm.results$estimate))
		# TT <- as.factor(as.character(rlm.results$Drug:rlm.results$Dose))
		# rlm.results$time <- relevel(rlm.results$time, ref="4h")

		aov1 <- aov(estimate~Drug+time, data=rlm.results)
		p.drug <- unlist(summary(aov1))[13]
		p.time <- unlist(summary(aov1))[14]

		# Graphe
		rlm.results <- rlm.results[order(rlm.results$time),]
		s.val <- rlm.results$error
		estim <- rlm.results$estimate

		col.palette <- c(rep("indianred1",2), rep("violetred",2))	#c(rep("plum",4), rep("mediumorchid1",4), rep("darkviolet",4))
		boxplot(rlm.results$estimate~rlm.results$Drug:rlm.results$time,  border=NA,
			ylim=range(min(estim)-1, max(estim)+1), names=rep(c("Placebo",Drug.test[D]), 2), 
			main=unique(rlm.results$Comp), sub = paste(Drug.test[D], "treatment"),
			cex.main = 2, cex.sub = 2, cex.axis = 1.5, par(bg = "ivory"))

		for (i in 1:nrow(rlm.results)){
			points(i, estim[i], pch=18, cex=4, col=col.palette[i])
			segments(x0=i, x1=i, y0=estim[i]-s.val[i]*2, y1=estim[i]+s.val[i]*2, lwd=3, col=col.palette[i])
			segments(x0=i-0.1, x1=i+0.1, y0=estim[i]-s.val[i]*2, y1=estim[i]-s.val[i]*2, lwd=3, col=col.palette[i])
			segments(x0=i-0.1, x1=i+0.1, y0=estim[i]+s.val[i]*2, y1=estim[i]+s.val[i]*2, lwd=3, col=col.palette[i])
			}

		for (i in seq(1,nrow(rlm.results), by=2))
			lines(estim[i:(i+1)]~seq(i, (i+1)), lwd=3, col=col.palette[i])

		legend1 <- paste("Drug effect: p-value =", signif(p.drug,3))	
		legend2 <- paste("Time effect: p-value =", signif(p.time,3))	

		text(x = c(0.5, 0.5), y=c(max(estim)+0.9, max(estim)+0.7), offset = 0, pos = 4, 
				labels=c(legend1, legend2), font=4, cex=1.5, col="royalblue4")
				

		print(list(res =rlm.results, Summary=summary(aov1)))

		All.Res <- rbind(All.Res, rlm.results)
			
	}

	return(list(All.Res = All.Res, summary(aov2<-aov(estimate~Drug+time, data=All.Res))))
}
