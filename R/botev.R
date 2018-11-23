# Compute the optimal kernel bandwidth
botev <- function(x){
    n <- 512
    minimum <- min(x,na.rm=TRUE)
    maximum <- max(x,na.rm=TRUE)
    Range <- maximum - minimum
    MIN <- minimum-Range/10
    MAX <- maximum+Range/10
    R <- MAX-MIN
    dx <- R/n
    xmesh <- MIN+seq(0,R,dx)
    if (anyDuplicated(x)){
        N <- length(as.numeric(names(table(x))))
    } else {
        N <- length(x)
    }
    w <- graphics::hist(x,xmesh,plot=FALSE)
    initialdata <- (w$counts)/N
    initialdata <- initialdata/sum(initialdata)
    a <- dct1d(initialdata)
    II <- (1:(n-1))^2
    a2 <- (a[2:n]/2)^2
    tstar <- tryCatch(stats::uniroot(fixedpoint,c(0,.1),
                      N=N,II=II,a2=a2,tol=10^(-22))$root,
                      error=function(e).28*N^(-2/5))
    bandwidth <- sqrt(tstar)*R
    return(bandwidth)
}

# Botev:
# computes the discrete cosine transform of the column vector data
dct1d <- function(dat){
    n <- length(dat)
    # Compute weights to multiply DFT coefficients
    weight <- c(1,2*exp(-1i*(1:(n-1))*pi/(2*n)))
    # Re-order the elements of the columns of x
    dat <- c(dat[seq(1,n-1,2)], dat[seq(n,2,-2)])
    # Multiply FFT by weights:
    out <- Re(weight*stats::fft(dat))
    return(out)
}

# Botev:
# this implements the function t-zeta*gamma^[l](t)
fixedpoint <-  function(tt,N,II,a2){
    l <- 7
    f <- 2*(pi^(2*l))*sum((II^l)*a2*exp(-II*(pi^2)*tt))
    for (s in (l-1):2){
        K0 <- prod(seq(1,2*s-1,2))/sqrt(2*pi)
        const <- (1+(1/2)^(s+1/2))/3
        TT <- (2*const*K0/N/f)^(2/(3+2*s))
        f <- 2*pi^(2*s)*sum(II^s*a2*exp(-II*pi^2*TT))
    }
    out <- tt-(2*N*sqrt(pi)*f)^(-2/5)
    return(out)
}

# returns the median bandwidth for plotting several KDEs together
commonbandwidth <- function(x,log=FALSE){
    dnames <- names(x)
    n <- length(dnames)
    bw <- rep(0,n)
    for (i in 1:n){
        if (log) bw[i] <- botev(log(x[[i]]))
        else bw[i] <- botev(x[[i]])
    }
    return(stats::median(bw))
}
