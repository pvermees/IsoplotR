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
    McL <- mclean(tt=exp(fit$par['t']),d=x$d)
    if (x$d$U48$option==2){
        lims$U48i <- c('ll'=log(McL$U48i)-2,'ul'=min(log(McL$U48i)+1,log(20)))
        out <- searchlimithelper(lims,'U48','U48i',init=fit$par,method=method,
                                 x=x,anchor=anchor,type=type,model=model)
    }
    if (x$d$ThU$option==2){
        lims$ThUi <- c('ll'=log(McL$ThUi)-2,'ul'=min(log(McL$ThUi)+1,log(20)))
        out <- searchlimithelper(lims,'ThU','ThUi',init=fit$par,method=method,
                                 x=x,anchor=anchor,type=type,model=model)
    }
    if ('t'%in%pnames){
        st <- sqrt(fit$cov['t','t'])
        out['ll','t'] <- out['ll','t'] - 2*st
        out['ul','t'] <- out['ul','t'] + 2*st
    }
    out[,colSums(is.na(out))==0]
}

bayeslud <- function(fit,x,anchor=0,type='joint',model=1,debug=FALSE){
    lims <- getsearchlimits(fit=fit,x=x,anchor=anchor,type=type,model=model)
    pnames <- colnames(lims)
    np <- length(pnames)
    parlist <- list()
    if ('t'%in%pnames){
        parlist$t <- seq(from=lims['ll','t'],to=lims['ul','t'],length.out=20)
    }
    if ('a0'%in%pnames){
        parlist$a0 <- seq(from=lims['ll','a0'],to=lims['ul','a0'],length.out=5)
    }
    if ('b0'%in%pnames){
        parlist$b0 <- seq(from=lims['ll','b0'],to=lims['ul','b0'],length.out=5)
    }
    if ('U48i'%in%pnames){
        parlist$U48i <- seq(from=lims['ll','U48i'],
                            to=lims['ul','U48i'],length.out=20)
    }
    if ('ThUi'%in%pnames){
        parlist$ThUi <- seq(from=lims['ll','ThUi'],
                            to=lims['ul','ThUi'],length.out=20)
    }
    ni <- lapply(parlist,'length')
    pargrid <- expand.grid(parlist)
    LLgrid <- array(NA,dim=unlist(lapply(parlist,'length')),dimnames=parlist)
    for (i in 1:nrow(pargrid)){
        p <- as.numeric(pargrid[i,])
        names(p) <- pnames
        LLgrid[i] <- LL.ludwig(p,x=x,model=model,
                               anchor=anchor,type=type,debug=FALSE)
    }
    Lgrid <- exp(min(LLgrid)-LLgrid) # scale to avoid machine precision issues
    out <- list()
    for (pname in pnames){
        out[[pname]] <- cbind(exp(parlist[[pname]]),
                              log(colSums(apply(Lgrid,pname,'+'))))
    }
    if (debug){
        plot(out$t,type='b')
    }
    out
}
