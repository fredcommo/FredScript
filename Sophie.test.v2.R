# Function
Sophie.test.v2 <- function(data, Test, Ref, Drug.test=c("SU","ZD")){

	sdata <- split(data, data$Replic)
	n <- length(Drug.test)
	par(mfrow=c(1,n))
	All.Res <- c()

	R <- length(sdata)
	
	for (r in 1:R){
	data <- sdata[[r]]
		for(D in 1:n){
			index.D <- which(data$Drug== Drug.test[D] | data$Drug== "Placebo")	# Valeurs pour ZD & Placebo
			sub_D <- data[index.D,]	#diff2[index.ZD,]
			test = Test
			ref = Ref
			p.prot = as.character(ABlist$Antibody[test])	
			prot = as.character(ABlist$Antibody[ref])					
			index.prot <- which(colnames(sub_D)==prot)		# which(fileInfo$Ab==prot)
			index.p.prot <- which(colnames(sub_D)==p.prot)		# which(fileInfo$Ab==prot)


			rlm.results <- c()
				for (t in c(2, 1, 3)){
					time = levels(data$Time)[t]	
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
						for (j in 1:9) tmp <- rbind.data.frame(tmp, D_time[i,1:5])	
				
						tmp <- cbind.data.frame(tmp, delta)
						ratio <- rbind.data.frame(ratio, tmp)	# tous les delta phospho vs. total, pour toute concentration, au temps T, inclus ZD & Placebo 
						}

					drug <- as.factor(as.character(ratio$Drug))
					dose <- as.factor(as.character(ratio$Dose))
					dd <- as.factor(as.character(drug:dose))
					dil <- ratio$Dilution

					# Estimateur de la valeur moyenne par rlm

					for (d in 1:nlevels(drug)){
					d.dose <- c("0,1", "0,5", "1")
					if (levels(drug)[d]=="Placebo") d.dose <- "0"
			
						for (dos in 1:length(d.dose)){
							test <- ratio$delta[which(drug==levels(drug)[d] & dose==d.dose[dos])]
							dil.test <- log2(dil[which(drug==levels(drug)[d] & dose==d.dose[dos])])

							rob <- lm.rob(dil.test, test) # ; rlmB <- rlm(testB~dilB)

							wi <- rob$weights
							s2 <- sum(wi*(predict(rob)-test)^2)/sum(wi)
							s <- sqrt(s2*(1+1/sum(wi)))

							newx = (min(dil.test) + max(dil.test))/2 #; newxB <- (min(dilB) + max(dilB))/2
							estim <-predict(rob, newdata=list(x=newx))
	
							rlm.results<-rbind(rlm.results, c(time, levels(drug)[d], d.dose[dos], names(sdata)[r], estim, s))
							}
					}
				}

			rlm.results <- as.data.frame(rlm.results)
			colnames(rlm.results) <- c("Time", "Drug", "Dose", "Replic", "estimate", "error")
			rlm.results$estimate <- as.numeric(as.vector(rlm.results$estimate))
			rlm.results$error <- as.numeric(as.vector(rlm.results$error))
			rlm.results <- cbind.data.frame(Comp=paste(p.prot," vs. ",prot, sep=""), rlm.results)

			# Analyse par Anova
			#estimate <- as.numeric(as.vector(rlm.results$estimate))
			# TT <- as.factor(as.character(rlm.results$Drug:rlm.results$Dose))
			rlm.results$Time <- relevel(rlm.results$Time, ref="4h")

			aov1 <- aov(estimate ~ Dose + Time, data=rlm.results)
			p.dose <- unlist(summary(aov1))[13]
			p.time <- unlist(summary(aov1))[14]

			# Graphe
			rlm.results <- rlm.results[order(rlm.results$Time, rlm.results$Dose),]
			s.val <- rlm.results$error
			estim <- rlm.results$estimate

			col.palette <- c(rep("indianred1",4), rep("slateblue3",4), rep("violetred",4))	#c(rep("plum",4), rep("mediumorchid1",4), rep("darkviolet",4))
			boxplot(rlm.results$estimate~rlm.results$Dose:rlm.results$Time,  border=NA,
				ylim=range(min(estim)-0.5, max(estim)+1),
				names=rlm.results$Dose, main=unique(rlm.results$Comp), 
				cex.main = 2, cex.sub = 2, cex.axis = 1.5, par(bg = "ivory"))

			for (i in 1:nrow(rlm.results)){
				points(i, estim[i], pch=18, cex=4, col=col.palette[i])
				segments(x0=i, x1=i, y0=estim[i]-s.val[i], y1=estim[i]+s.val[i], lwd=3, col=col.palette[i])
				segments(x0=i-0.1, x1=i+0.1, y0=estim[i]-s.val[i], y1=estim[i]-s.val[i], lwd=3, col=col.palette[i])
				segments(x0=i-0.1, x1=i+0.1, y0=estim[i]+s.val[i], y1=estim[i]+s.val[i], lwd=3, col=col.palette[i])
				}

			for (i in seq(1,nrow(rlm.results), by=4))
				lines(estim[i:(i+3)]~seq(i, (i+3)), lwd=3, col=col.palette[i])

			legend1 <- paste("Dose effect: p-value =", signif(p.dose,3))	
			legend2 <- paste("Time effect: p-value =", signif(p.time,3))	

			text(x = c(0.5, 0.5), y=c(max(estim)+0.9, max(estim)+0.7), offset = 0, pos = 4, 
					labels=c(legend1, legend2), font=4, cex=1.5, col="royalblue4")
				

			print(list(res =rlm.results, Summary=summary(aov1)))

			All.Res <- rbind(All.Res, rlm.results)
			
		}
	}
	if(length(Drug.test)==1) aov2<-aov(estimate~Dose*Time+Replic, data=All.Res)
	else aov2<-aov(estimate~Drug+Dose+Time, data=All.Res)

	return(list(All.Res = All.Res, summary(aov2)))
}
