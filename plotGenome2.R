plotGenome2 <- function (aCGH.obj, samples = 1:num.samples(aCGH.obj), naut = 22, 
    Y = TRUE, X = TRUE, data = log2.ratios(aCGH.obj), chrominfo = human.chrom.info.Jul03, 
    yScale = c(-2, 2), samplenames = sample.names(aCGH.obj), ylb = "Log2Ratio") 
{
    datainfo <- clones.info(aCGH.obj)
    nchr <- naut
    if (X) 
        nchr <- nchr + 1
    if (Y) 
        nchr <- nchr + 1
    nsamples <- length(samplenames)
    ord <- order(datainfo$Chrom, datainfo$kb)
    chrom <- datainfo$Chrom[ord]
    kb <- datainfo$kb[ord]
    data <- data[ord, ]
    ind.unmap <- which(is.na(chrom) | is.na(kb) | (chrom > (naut + 2)))
    if (length(ind.unmap) > 0) {
        chrom <- chrom[-ind.unmap]
        kb <- kb[-ind.unmap]
        data <- data[-ind.unmap, ]
    }
    data <- data[chrom <= nchr, ]
    kb <- kb[chrom <= nchr]
    chrom <- chrom[chrom <= nchr]
    chrominfo <- chrominfo[1:nchr, ]
    chrom.start <- c(0, cumsum(chrominfo$length))[1:nchr]
    chrom.centr <- chrom.start + chrominfo$centr
    chrom.mid <- chrom.start + chrominfo$length/2
    chrom.rat <- chrominfo$length/max(chrominfo$length)
    par(cex = 0.6, pch = 18, lab = c(1, 6, 7), cex.axis = 1.5, xaxs = "i")

    for (k in 1:length(samples)) {
        vec <- data[, samples[k]]
        name <- samplenames[samples[k]]
        clone.genomepos <- rep(0, length(kb))
        for (i in 1:nrow(chrominfo)) clone.genomepos[chrom == i] <- kb[chrom == i] + chrom.start[i]
        y.ranges <- sapply(1:nrow(chrominfo), function(i) range(vec[chrom == i], yScale, na.rm = TRUE))
        ylim <- range(y.ranges)

        plot(clone.genomepos/1000, vec, ylim = ylim, xlab = "", 
            ylab = "", xlim = c(min(clone.genomepos[clone.genomepos > 0], 
		na.rm = TRUE)/1000, clone.genomepos[sum(clone.genomepos > 0)]/1000), 
		col=ifelse(vec<0, "firebrick2", ifelse(vec>0, "green3", "grey10")), xaxt = "n")

        axis(side = 1, at = clone.genomepos[1]/1000, label = "", tick = FALSE)
        title(main = paste(samples[k], " ", name), ylab = ylb, xlab = "", cex.lab = 1.5, cex.main = 2)
        for (i in seq(1, naut, b = 2)) mtext(paste("", i), side = 1, at = chrom.mid[i]/1000, line = 0.3, col = "red")
        for (i in seq(2, naut, b = 2)) mtext(paste("", i), side = 3, at = chrom.mid[i]/1000, line = 0.3, col = "red")
        if (X) mtext("X", side = 1, at = chrom.mid[naut + 1]/1000, line = 0.3, col = "red")
        if (Y) mtext("Y", side = 3, at = chrom.mid[naut + 2]/1000, line = 0.3, col = "red")

        abline(v = c(chrom.start/1000, (chrom.start[nrow(chrominfo)] + chrominfo$length[nrow(chrominfo)])/1000), lty = 1, col = "grey20")
        abline(h = seq(-1, 1, b = 0.5), lty = 3, col = "grey40")
        abline(v = (chrominfo$centromere + chrom.start)/1000, lty = 3, col = "royalblue1")
    }
    invisible(list(x = clone.genomepos/1000, ylim = range(y.ranges)))
}
