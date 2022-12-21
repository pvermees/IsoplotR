recursivelimitsearch <- function(pname,iname,ll,ul,LLmax,method,
                                 x=x,anchor=0,type=1,model=1,
                                 maxlevel=5,side='lower',debug=FALSE){
    newlim <- (ll+ul)/2
    X <- x
    X$d[[pname]] <- list(x=newlim,sx=0,option=1)
    init <- init.ludwig(x=X,model=model,anchor=anchor,type=type)
    fit <- stats::optim(init,fn=LL.ludwig,method=method,hessian=FALSE,
                        x=X,anchor=anchor,type=type,model=model)
    par <- fit$par
    par[iname] <- log(newlim)
    LL <- LL.ludwig(par,x=x,anchor=anchor,type=type,model=model)
    if (debug) print(c(LLmax,LL,ll,ul))
    out <- list(ta0b0=exp(fit$par))
    if (LL<(LLmax+100)){ # not far enough
        if (maxlevel<1){
            if (side=='lower'){
                out[[pname]] <- ll
            } else {
                out[[pname]] <- ul
            }
            return(out)
        } else {
            if (side=='lower'){
                ul <- newlim
            } else {
                ll <- newlim
            }
        }
    } else { # too far
        if (maxlevel<1){
            out[[pname]] <- newlim
            return(out)
        } else {
            if (side=='lower'){
                ll <- newlim
            } else {
                ul <- newlim
            }
        }        
    }
    recursivelimitsearch(pname=pname,iname=iname,ll=ll,ul=ul,LLmax=LLmax,
                         method=method,x=x,anchor=anchor,type=type,
                         model=model,maxlevel=maxlevel-1,side=side)
}
    
searchlimithelper <- function(lims,pname,iname,fit,method,x=x,
                              anchor=0,type=1,model=1,debug=FALSE){
    McL <- mclean(tt=exp(fit$par['t']),d=x$d)
    par <- fit$par
    par[pname] <- log(unname(McL$U48i))
    LLmax <- LL.ludwig(par,x=x,anchor=anchor,type=type,model=model)
    if (debug){
        U48i <- seq(from=0,to=9,length.out=20)
        LL <- U48i*0
        for (i in seq_along(U48i)){
            X <- x
            X$d$U48 <- list(x=U48i[i],sx=0,option=1)
            init <- init.ludwig(x=X,model=model,anchor=anchor,type=type)
            fit <- stats::optim(init,fn=LL.ludwig,method=method,hessian=FALSE,
                                x=X,anchor=anchor,type=type,model=model)
            par <- fit$par
            par['U48i'] <- log(U48i[i])
            LL[i] <- LL.ludwig(par,x=x,anchor=anchor,type=type,model=model)
        }
        plot(U48i,LL,type='b')
        lines(c(0,20),rep(LLmax,2))
        lines(rep(McL$U48i,2),range(LL))
    }
    llfit <- recursivelimitsearch(pname=pname,iname=iname,ll=0,ul=unname(McL$U48i),
                                  LLmax=LLmax,method=method,x=x,anchor=anchor,
                                  type=type,model=model,maxlevel=5,side='lower')
    ulfit <- recursivelimitsearch(pname=pname,iname=iname,ll=unname(McL$U48i),ul=20,
                                  LLmax=LLmax,method=method,x=x,anchor=anchor,
                                  type=type,model=model,maxlevel=5,side='upper')
    out <- lims
    out['ll',iname] <- llfit[[pname]]
    out['ul',iname] <- ulfit[[pname]]
    pnames <- colnames(lims)
    if ('t'%in%pnames){
        out['ll','t'] <- min(lims['ll','t'],llfit$ta0b0['t'],
                             ulfit$ta0b0['t'],na.rm=TRUE)
        out['ul','t'] <- max(lims['ul','t'],llfit$ta0b0['t'],
                             ulfit$ta0b0['t'],na.rm=TRUE)
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
        out <- searchlimithelper(lims,pname='U48',iname='U48i',fit=fit,
                                 method=method,x=x,anchor=anchor,
                                 type=type,model=model)
    }
    if (x$d$ThU$option==2){
        out <- searchlimithelper(lims,pname='ThU',iname='ThUi',fit=fit,
                                 method=method,x=x,anchor=anchor,
                                 type=type,model=model)
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
    ilist <- list()
    if ('t'%in%pnames){
        nt <- 25
        parlist$t <- seq(from=lims['ll','t'],to=lims['ul','t'],length.out=nt)
        ilist$t <- 1:nt
    }
    if ('a0'%in%pnames){
        na0 <- 7
        parlist$a0 <- seq(from=lims['ll','a0'],to=lims['ul','a0'],length.out=na0)
        ilist$a0 <- 1:na0
    }
    if ('b0'%in%pnames){
        nb0 <- 7
        parlist$b0 <- seq(from=lims['ll','b0'],to=lims['ul','b0'],length.out=nb0)
        ilist$b0 <- 1:nb0
    }
    if ('U48i'%in%pnames){
        n48i <- 25
        parlist$U48i <- seq(from=lims['ll','U48i'],
                            to=lims['ul','U48i'],length.out=n48i)
        ilist$U48i <- 1:n48i
    }
    if ('ThUi'%in%pnames){
        nThUi <- 25
        parlist$ThUi <- seq(from=lims['ll','ThUi'],
                            to=lims['ul','ThUi'],length.out=nThUi)
        ilist$ThUi <- 1:nThUi
    }
    ni <- lapply(parlist,'length')
    LLgrid <- array(NA,dim=unlist(lapply(parlist,'length')),dimnames=parlist)
    logparlist <- lapply(parlist,'log')
    logpargrid <- expand.grid(logparlist)
    igrid <- expand.grid(ilist)
    for (i in 1:nrow(logpargrid)){
        p <- as.numeric(logpargrid[i,])
        names(p) <- pnames
        LLgrid[i] <- -LL.ludwig(p,x=x,model=model,anchor=anchor,type=type)
    }
    out <- list()
    for (pname in pnames){
        x <- parlist[[pname]]
        y <- 0*x
        for (i in seq_along(x)){
            j <- which(igrid[,pname]%in%ilist[[pname]][i])
            y[i] <-log_sum_exp(u=LLgrid[j[1]],v=LLgrid[j[-1]])
        }
        out[[pname]] <- cbind(x,y)
    }
    if (debug){
        bayesplot <- function(dat,...){
            opt <- 2
            if (opt==1){
                xy <- bayespline(dat)
                #plot(xy$x,cumsum(exp(xy$y)),type='l',...)
                plot(xy$x,xy$y,type='l',...)
            } else if (opt==2){
                x <- dat[,1]
                y <- exp(dat[,2]-max(dat[,2]))
                plot(x,y,...,type='b')
            } else {
                plot(dat,...,type='b')
            }
        }
        op <- par(mfrow=c(2,2))
        if ('t'%in%pnames) bayesplot(out$t,xlab='t',ylab='F(t)')
        if ('a0'%in%pnames) bayesplot(out$a0,xlab='a0',ylab='F(a0)')
        if ('U48i'%in%pnames) bayesplot(out$U48i,xlab='U48i',ylab='F(U48i)')
        if ('ThUi'%in%pnames) bayesplot(out$ThUi,xlab='ThUi',ylab='F(ThUi)')
        par(op)
    }
    out
}

bayespline <- function(xy,n=100){
    good <- is.finite(xy[,2])
    spline(xy[good,1],xy[good,2],n=n)
}
