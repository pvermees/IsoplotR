compisofit <- function(dat){
    
    logcovmat <- function(X,sX,Y,sY,sXY){
        ns <- length(sX)
        E <- matrix(0,2*ns,2*ns)
        iX <- 1:ns
        iY <- ns+iX
        diag(E)[iX] <- sX^2
        diag(E)[iY] <- sY^2
        diag(E[iX,iY]) <- sXY
        diag(E[iY,iX]) <- sXY
        J <- matrix(0,2*ns,2*ns)
        diag(J)[iX] <- 1/X
        diag(J)[iY] <- 1/Y
        J %*% E %*% t(J)
    }
    
    l1xx0misfit <- function(l1xx0,lx0,ly0,X,Y,O,LL=TRUE){
        x0 <- exp(lx0)
        y0 <- exp(ly0)
        x <- x0*exp(1-exp(l1xx0))
        ns <- length(X)
        D <- rep(NA,ns)
        iX <- 1:ns
        iY <- ns+iX
        D[iX] <- log(X)-log(x)
        D[iY] <- log(Y)-ly0-l1xx0
        SS <- D %*% O %*% D
        if (LL) out <- SS/2
        else out <- SS
        out
    }
    
    misfit <- function(lxy0,X,Y,O,LL=TRUE){
        x0 <- exp(lxy0[1])
        il1xx0 <- log(abs(1-X/x0))
        l1xx0 <- optim(il1xx0,fn=l1xx0misfit,method='BFGS',
                       lx0=lxy0[1],ly0=lxy0[2],X=X,Y=Y,O=O)$par
        l1xx0misfit(l1xx0,lx0=lxy0[1],ly0=lxy0[2],X=X,Y=Y,O=O,LL=LL)
    }
    
    X <- dat[,1]
    sX <- dat[,2]
    Y <- dat[,3]
    sY <- dat[,4]
    rXY <- dat[,5]
    sXY <- rXY*sX*sY
    E <- logcovmat(X,sX,Y,sY,sXY)
    O <- solve(E)
    
    ab <- abs(lm(Y ~ X)$coefficients)
    init <- c(log(ab[1]/ab[2]),log(ab[1]))
    fit <- optim(init,fn=misfit,X=X,Y=Y,O=O,hessian=TRUE)
    out <- list()
    out$logpar <- fit$par
    out$logcov <- solve(fit$hessian)
    out$LL <- fit$value
    out$SS <- misfit(fit$par,X=X,Y=Y,O=O,LL=FALSE)
    out
}
