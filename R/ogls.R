LL.ogls <- function(ab,dat,omega){
    a <- ab[1]
    b <- ab[2]
    ns <- nrow(dat)/2
    iX <- seq(from=1,to=2*ns-1,length.out=ns)
    iY <- seq(from=2,to=2*ns,length.out=ns)
    X <- dat[iX,1,drop=FALSE]
    Y <- dat[iY,1,drop=FALSE]
    rY <- Y - a - b * X
    AA <- omega[iX,iX]+2*b*omega[iX,iY]+omega[iY,iY]*b^2
    BB <- (omega[iX,iY]+b*omega[iY,iY])%*%rY
    rx <- - solve(AA,BB)
    fitted.x <- X - rx
    ry <- dat[iY,1] - a - b*fitted.x
    v <- matrix(c(rx,ry),nrow=1,ncol=2*ns)
    (v %*% omega %*% t(v))
}
ogls <- function(dat){
    out <- list()
    ydat <- data2york.ogls(dat)
    yfit <- york(ydat)
    init <- c(yfit$a[1],yfit$b[1])
    omega <- solve(dat[,-1])
    fit <- optim(init,LL.ogls,dat=dat,omega=omega,hessian=TRUE)
    covmat <- solve(fit$hessian)
    out$a <- c(fit$par[1],sqrt(covmat[1,1]))
    out$b <- c(fit$par[2],sqrt(covmat[2,2]))
    out$cov.ab <- covmat[1,2]
    out
}
