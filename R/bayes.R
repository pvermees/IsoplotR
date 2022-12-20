searchlimithelper <- function(lims,pname,cname,init,m=0,M=20,method,x=x,
                              anchor=0,type=1,model=1){
    out <- lims
    X <- x
    X$d[[pname]] <- list(x=lims['ll',cname],sx=0,option=1)
    llfit <- stats::optim(init,fn=LL.ludwig,method=method,hessian=FALSE,
                          x=X,anchor=anchor,type=type,model=model)$par
    X$d[[pname]] <- list(x=lims['ul',cname],sx=0,option=1)
    ulfit <- stats::optim(init,fn=LL.ludwig,method=method,hessian=FALSE,
                          x=X,anchor=anchor,type=type,model=model)$par
    pnames <- colnames(lims)
    if ('t'%in%pnames){
        out['ll','t'] <- min(lims['ll','t'],llfit['t'],ulfit['t'],na.rm=TRUE)
        out['ul','t'] <- max(lims['ul','t'],llfit['t'],ulfit['t'],na.rm=TRUE)
    }
    if ('a0'%in%pnames){
        out['ll','a0'] <- min(lims['ll','a0'],llfit['a0'],ulfit['a0'],na.rm=TRUE)
        out['ul','a0'] <- max(lims['ul','a0'],llfit['a0'],ulfit['a0'],na.rm=TRUE)
    }
    if ('b0'%in%pnames){
        out['ll','b0'] <- min(lims['ll','b0'],llfit['b0'],ulfit['b0'],na.rm=TRUE)
        out['ul','b0'] <- max(lims['ul','b0'],llfit['b0'],ulfit['b0'],na.rm=TRUE)
    }
    out
}

getsearchlimits <- function(fit,x,anchor=0,type='joint',model=1){
    pnames <- names(fit$par)
    np <- length(pnames)
    if (np>1) method <- "Nelder-Mead"
    else method <- "BFGS"
    lims <- matrix(NA,nrow=2,ncol=np)
    rownames(lims) <- c('ll','ul')
    colnames(lims) <- pnames
    lims <- as.data.frame(lims)
    if (x$d$U48$option==2){
        lims$U48i <- c('ll'=0,'ul'=20)
        out <- searchlimithelper(lims,'U48','U48i',init=fit$par,method=method,
                                 x=x,anchor=anchor,type=type,model=model)
    }
    if (x$d$ThU$option==2){
        out <- searchlimithelper(lims,'ThU','ThUi',init=fit$par,method=method,
                                 x=x,anchor=anchor,type=type,model=model)
    }
    out
}

bayeslud <- function(fit,x,anchor=0,type='joint',model=1){
    lims <- getsearchlimits(fit=fit,x=x,anchor=anchor,type=type,model=model)
}
