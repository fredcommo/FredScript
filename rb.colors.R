rb.colors<-function (n, alpha = 1) 
{
    if ((n <- as.integer(n[1])) > 0) {
        even.n <- n%%2 == 0
        k <- n%/%2
        l1 <- k + 1 - even.n
        l2 <- n - k + even.n
        c(if (l1 > 0) hsv(h =11/20, s = seq.int(0.95, ifelse(even.n, 
            0.95/k, 0), length.out = l1), v = 0.75, alpha = alpha), 
            if (l2 > 1) hsv(h = 20/20, s = seq.int(0, 0.9, length.out = l2)[-1], 
                v = 0.75, alpha = alpha))
    }
    else character(0)
}

# Gamme couleurs : image(z=matrix(seq(1:100),100,1),col=rb.colors(100,1))

# Gamme couleurs : 
# x = seq(0, 10, length=10000)
# y = seq(0,1,by=0.1)
# image(x=x, z=matrix(seq(0, 10, length=9999),length(x)-1,1), xlim=range(0, 10),col=rainbow(length(delta), start = 0.7, end = 0.3))

