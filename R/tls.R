# total least squares using PCA, with jackknife error estimation
# dat = data matrix whose first column is independent variable
tls <- function(dat){
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
    out
}

tlspar <- function(dat){
    pc <- stats::prcomp(dat)
    np <- ncol(dat)-1
    out <- vector()
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
