compisofit <- function(x,type=1){
    
    logcovmat <- function(X,sX,Y,sY,sXY){
        ns <- length(sX)
        E <- matrix(0,2*ns,2*ns)
        iX <- seq(from=1,to=2*ns-1,by=2)
        iY <- iX + 1
        diag(E)[iX] <- sX^2
        diag(E)[iY] <- sY^2
        diag(E[iX,iY]) <- sXY
        diag(E[iY,iX]) <- sXY
        J <- matrix(0,2*ns,2*ns)
        diag(J)[iX] <- 1/X
        diag(J)[iY] <- 1/Y
        J %*% E %*% t(J)
    }
    
    misfit <- function(lx,lx0,ly0,X,Y,O){
        x <- exp(lx)
        x0 <- exp(lx0)
        y0 <- exp(ly0)
        ns <- length(X)
        D <- rep(NA,ns)
        iX <- seq(from=1,to=2*ns-1,by=2)
        iY <- iX + 1
        D[iX] <- log(X)-log(x)
        D[iY] <- log(Y)-log(y0)-log(1-x/x0)
        D %*% O %*% D
    }
    
    xy0misfit <- function(lxy0,X,Y,O){
        lx <- optim(log(X),misfit,lx0=lxy0[1],ly0=lxy0[2],
                    X=X,Y=Y,O=O,method='BFGS')$par
        misfit(lx,lx0=lxy0[1],ly0=lxy0[2],X=X,Y=Y,O=O)
    }
    
    X <- dat[,1]
    sX <- dat[,2]
    Y <- dat[,3]
    sY <- dat[,4]
    rXY <- dat[,5]
    sXY <- rXY*sX*sY
    E <- logcovmat(X,sX,Y,sY,sXY)
    O <- solve(E)
    
    ab <- lm(Y ~ X)$coefficients
    init <- c(log(-ab[1]/ab[2]),log(ab[1]))
    lxy0fit <- optim(init,xy0misfit,X=X,Y=Y,O=O,method='BFGS')$par
    lxfit <- optim(log(X),misfit,lx0=lxy0fit[1],ly0=lxy0fit[2],
                   X=X,Y=Y,O=O,method='BFGS')$par
    
    out <- list()
    out$x <- exp(lxfit)
    out$x0 <- exp(lxy0fit)[1]
    out$y0 <- exp(lxy0fit)[2]
    out
}
