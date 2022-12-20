searchlimithelper <- function(lims,pname,cname,init,m=0,M=20,method,x=x,
                              anchor=0,type=1,model=1){
    out <- lims
    X <- x
    X$d[[pname]] <- list(x=lims['ll',cname],sx=0,option=1)
    llfit <- exp(stats::optim(init,fn=LL.ludwig,method=method,hessian=FALSE,
                              x=X,anchor=anchor,type=type,model=model)$par)
    X$d[[pname]] <- list(x=lims['ul',cname],sx=0,option=1)
    ulfit <- exp(stats::optim(init,fn=LL.ludwig,method=method,hessian=FALSE,
                              x=X,anchor=anchor,type=type,model=model)$par)
    pnames <- colnames(lims)
    if ('t'%in%pnames){
        out['ll','t'] <- min(lims['ll','t'],llfit['t'],ulfit['t'],na.rm=TRUE)
        out['ul','t'] <- max(lims['ul','t'],llfit['t'],ulfit['t'],na.rm=TRUE)
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
        lims$U48i <- c('ll'=0,'ul'=min(max(1,2*McL$U48i),20))
        out <- searchlimithelper(lims,'U48','U48i',init=fit$par,method=method,
                                 x=x,anchor=anchor,type=type,model=model)
    }
    if (x$d$ThU$option==2){
        lims$ThUi <- c('ll'=0,'ul'=min(max(1,2*McL$ThUi),20))
        out <- searchlimithelper(lims,'ThU','ThUi',init=fit$par,method=method,
                                 x=x,anchor=anchor,type=type,model=model)
    }
    if ('t'%in%pnames){
        st <- sqrt(fit$cov['t','t'])
        dt <- log(diff(out[,'t']))
        out['ll','t'] <- exp(log(out['ll','t']) - max(2*st,dt/10))
        out['ul','t'] <- exp(log(out['ul','t']) + max(2*st,dt/10))
    }
    if ('a0'%in%pnames){
        sa0 <- sqrt(fit$cov['a0','a0'])
        out['ll','a0'] <- exp(fit$par['a0'] - 5*sa0)
        out['ul','a0'] <- exp(fit$par['a0'] + 5*sa0)
    }
    if ('b0'%in%pnames){
        sb0 <- sqrt(fit$cov['b0','b0'])
        out['ll','b0'] <- exp(fit$par['b0'] - 5*sb0)
        out['ul','b0'] <- exp(fit$par['b0'] + 5*sb0)
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
        parlist$a0 <- seq(from=lims['ll','a0'],to=lims['ul','a0'],length.out=7)
    }
    if ('b0'%in%pnames){
        parlist$b0 <- seq(from=lims['ll','b0'],to=lims['ul','b0'],length.out=7)
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
    LLgrid <- array(NA,dim=unlist(lapply(parlist,'length')),dimnames=parlist)
    logparlist <- lapply(parlist,'log')
    logpargrid <- expand.grid(logparlist)
    for (i in 1:nrow(logpargrid)){
        p <- as.numeric(logpargrid[i,])
        names(p) <- pnames
        LLgrid[i] <- LL.ludwig(p,x=x,model=model,
                               anchor=anchor,type=type,debug=FALSE)
    }
    Lgrid <- exp(min(LLgrid)-LLgrid) # scale to avoid machine precision issues
    out <- list()
    for (pname in pnames){
        out[[pname]] <- cbind(parlist[[pname]],
                              log(colSums(apply(Lgrid,pname,'+'))))
    }
    if (debug){
        bayesplot <- function(xy,...){
            plot(xy$x,cumsum(exp(xy$y)),type='l',...)
        }
        op <- par(mfrow=c(2,2))
        if ('t'%in%pnames) bayesplot(bayespline(out$t),
                                     xlab='t',ylab='F(t)')
        if ('a0'%in%pnames) bayesplot(bayespline(out$a0),
                                      xlab='a0',ylab='F(a0)')
        if ('U48i'%in%pnames) bayesplot(bayespline(out$U48i),
                                        xlab='U48i',ylab='F(U48i)')
        if ('ThUi'%in%pnames) bayesplot(bayespline(out$ThUi),
                                        xlab='ThUi',ylab='F(ThUi)')
        par(op)
    }
    out
}

bayespline <- function(xy,n=100){
    good <- is.finite(xy[,2])
    spline(xy[good,1],xy[good,2],n=n)
}
