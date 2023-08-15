LL.ogls <- function(ab,dat){
    get.fitted.x <- function(dat,a,b){
        ns <- length(dat$x)/2
        ix <- 1:ns
        iy <- (ns+1):(2*ns)
        X <- matrix(dat$x[ix],ns,1)
        Y <- matrix(dat$x[iy],ns,1)
        ry <- Y - a - b * X
        O <- dat$omega
        rx <- - solve(O[ix,ix]+2*b*O[ix,iy]+O[iy,iy]*b^2,
        (O[ix,iy]+b*O[iy,iy])%*%ry)
        out <- X - rx
        out
    }
    ns <- length(dat$x)/2
    a <- ab[1]
    b <- ab[2]
    fitted.x <- get.fitted.x(dat,a=a,b=b)
    rx <- dat$x[1:ns] - fitted.x[1:ns]
    ry <- dat$x[(ns+1):(2*ns)] - a - b*fitted.x[1:ns]
    v <- matrix(c(rx,ry),nrow=1,ncol=2*ns)
    (v %*% dat$omega %*% t(v))
}
ogls <- function(dat){
    out <- list()
    yorkdat <- data2york(dat)
    yfit <- york(yorkdat)
    out$init <- c(yfit$a[1],yfit$b[1])
    dat$omega <- solve(dat$covmat)
    fit <- optim(out$init,LL.ogls,dat=dat,hessian=TRUE)
    covmat <- solve(fit$hessian)
    out$a <- c(fit$par[1],sqrt(covmat[1,1]))
    out$b <- c(fit$par[2],sqrt(covmat[2,2]))
    out$cov.ab <- covmat[1,2]
    out
}
