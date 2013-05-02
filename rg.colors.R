rg.colors<-function (n, alpha = 1, Gamma = 5) 
{
    if ((n <- as.integer(n[1])) > 0) {
        even.n <- n%%2 == 0
        k <- n%/%2
        l1 <- k + 1 - even.n
        l2 <- n - k + even.n
        c(if (l1 > 0) hsv(h =8/20, s = seq.int(1, ifelse(even.n, 
            1/k, 0), length.out = l1), gamma = Gamma, v = 0.75, alpha = alpha), 
            if (l2 > 1) hsv(h = 20/20, s = seq.int(0, 1, length.out = l2)[-1], 
                gamma = Gamma, v = 0.75, alpha = alpha))
    }
    else character(0)
}


# Gamme couleurs : image(z=matrix(seq(1:100),100,1),col=rg.colors(100,1))
