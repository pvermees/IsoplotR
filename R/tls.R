# total least squares using PCA, with jackknife error estimation
# dat = data matrix whose first column is independent variable
tls <- function(dat,anchor=0){
    if (anchor[1]>0){
        out <- anchored.deming(dat,anchor=anchor)
    } else {
        out <- list()
        out$par <- tlspar(dat)
        ns <- nrow(dat)
        nv <- ncol(dat)
        np <- 2*(nv-1)
        jack <- matrix(0,ns,np)
        for (i in 1:ns){
            jack[i,] <- tlspar(dat[-i,])
        }
        avgjack <- colMeans(jack)
        jackdiff <- sweep(jack,2,avgjack,'-')
        SS <- matrix(0,np,np)
        for (i in 1:ns){
            SS <- SS + jackdiff[i,] %*% t(jackdiff[i,])
        }
        out$cov <- SS*(ns-1)/ns
        rownames(out$cov) <- colnames(out$cov) <- names(out$par)
    }
    out
}

tlspar <- function(dat){
    out <- vector()
    pc <- stats::prcomp(dat)
    np <- ncol(dat)-1
    for (i in 1:np){
        if (i%%2) pnames <- letters[2*i-c(1,0)]
        else pnames <- LETTERS[2*i-c(3,2)]
        slope <- pc$rotation[i+1,'PC1']/pc$rotation[1,'PC1']
        intercept <- pc$center[i+1] - slope*pc$center[1]
        out[pnames[1]] <- intercept
        out[pnames[2]] <- slope
    }
    out
}

anchored.deming <- function(dat,anchor){
    out <- list()
    out$par <- c(NA,NA)
    out$cov <- matrix(0,2,2)
    if (anchor[1]%in%c(1,'intercept','a') && length(anchor)>1){
        a <- anchor[2]
        init <- unname(stats::lm(I(dat[,2]-a) ~ 0 + dat[,1])$coefficients)
        interval <- sort(init*c(1/5,5))
        fit <- stats::optimise(deming.misfit.b,interval=interval,a=a,dat=dat)
        out$par <- c('a'=a,'b'=fit$minimum)
        H <- stats::optimHess(fit$minimum,deming.misfit.b,a=a,dat=dat)
        ve <- stats::var(deming_residuals(ab=out$par,dat=dat))
        out$cov[2,2] <- inverthess(H)*ve
    } else if (anchor[1]%in%c(2,'slope','b') && length(anchor)>1){
        b <- anchor[2]
        init <- unname(mean(dat[,2]-b*dat[,1]))
        interval <- sort(init*c(1/5,5))
        fit <- stats::optimise(deming.misfit.a,interval=interval,b=b,dat=dat)
        out$par <- c('a'=fit$minimum,'b'=b)
        H <- stats::optimHess(fit$minimum,deming.misfit.a,b=b,dat=dat)
        ve <- stats::var(deming_residuals(ab=out$par,dat=dat))
        out$cov[1,1] <- inverthess(H)*ve
    } else {
        stop("Invalid anchor.")
    }
    out
}

deming.misfit.a <- function(a,b,dat){
    deming.misfit.ab(c(a,b),dat)
}
deming.misfit.b <- function(b,a,dat){
    deming.misfit.ab(c(a,b),dat)
}
deming.misfit.ab <- function(ab,dat){
    if (ncol(dat)>2) stop("Anchored TLS regression currently only works in 2D")
    e <-deming_residuals(ab,dat) 
    sum(e^2)
}
deming_residuals <- function(ab,dat){
    a <- ab[1]
    b <- ab[2]
    x <- dat[,1]
    y <- dat[,2]    
    (y-a-b*x)/sqrt(1+b^2)    
}
