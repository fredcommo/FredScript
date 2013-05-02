
AddSegments <- function(seg.cna.obj, cgh, cutGainLoss = 10){

	# Ajoute les valeurs de segmentation pour chaque sonde/position

data <- cgh$cgh

seg.start <- seg.cna.obj$output$loc.start
seg.end <- seg.cna.obj$output$loc.end			# !!!! vérifier seg-end
seg.mean <- seg.cna.obj$output$seg.mean
seg.len <- seg.cna.obj$output$num.mark
cum.seg.len <- cumsum(seg.len)
s <- data$ChrStart[cum.seg.len]
e <- data$ChrEnd[cum.seg.len]
if(!is.null(e)) seg.end <- seg.end + (e-s)

# Proportion d'aberrations
cutoff <- log2(1 + cutGainLoss/100)
seg.delta <- seg.end - seg.start
Prop <- sum(seg.delta[which(abs(seg.mean)>=cutoff)])/sum(seg.delta)*100
Prop <- round(Prop, 2)

seg.mean <- seg.cna.obj$output$seg.mean
seg.val <- rep(0, nrow(data))
for(i in 1:length(seg.mean)){
	index <- which(data$GenomicPos>= seg.start[i] & data$GenomicPos<=seg.end[i])
	seg.val[index] <- seg.mean[i]
	}
cgh$cgh <- cbind.data.frame(data, Segm = seg.val)
cgh$Prop <- Prop
return(cgh)
}
