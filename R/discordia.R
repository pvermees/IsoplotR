concordia.age <- function(x,wetherill=TRUE){
    X <- UPb.preprocess(x,wetherill)
    fit <- optim(c(1,1), LL.concordia, x=X, method="BFGS", hessian=TRUE)
    out <- list()
    out$x <- fit$par
    out$cov <- solve(fit$hessian)
    names(out$x) <- names(X)
    colnames(out$cov) <- names(X)
    rownames(out$cov) <- names(X)
    out
}

UPb.preprocess <- function(x,wetherill){
    if (wetherill) selection <- c('Pb207U235','Pb206U238')
    else selection <- c('U238Pb206','Pb207Pb206')
    out <- list()
    for (i in 1:nrow(x$x)){
        X <- x$x[i,selection]
        covmat <- get.covmat.UPb(x,i)[selection,selection]
        out[[i]] <- list(x=X, cov=covmat)
    }
    out    
}

LL.concordia <- function(mu,x){
    LL <- 0
    for (i in 1:length(x)){
        X <- matrix(x[[i]]$x-mu,1,2)
        covmat <- x[[i]]$cov
        ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@37@"]]));##:ess-bp-end:##
        LL <- LL - log(2*pi) - 0.5*determinant(covmat,logarithmic=TRUE)$modulus -
                               0.5* X %*% solve(covmat) %*% t(X)
    }
    -LL
}
