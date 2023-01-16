initial2time <- function(x,anames,avalues,anchor=0,type='joint',
                         model=1,alim=c(0,20,1e-5),debug=FALSE){
    if (debug) browser()
    X <- x
    for (i in seq_along(anames)){
        X$d[[anames[i]]] <- list(x=avalues[i],sx=0,option=1)
    }
    init <- init.ludwig(x=X,model=model,anchor=anchor,type=type)
    fit <- stats::optim(init$par,fn=LL.ludwig,method='L-BFGS-B',
                        lower=init$lower,upper=init$upper,hessian=FALSE,
                        x=x,X=X,anchor=anchor,type=type,model=model,alim=alim)
    lt <- ifelse('t'%in%names(fit$par),fit$par['t'],anchor[2])
    McL <- mclean(tt=exp(lt),d=X$d)
    out <- list()
    out$par <- fit$par
    for (aname in anames){
        iname <- paste0(aname,'i')
        out$par[iname] <- McL[[iname]]
    }
    out$LL <- fit$value
    out
}
recursivelimitsearch_a <- function(aname,ll,ul,LLmax,x=x,anchor=0,type=1,
                                   model=1,LLbuffer=10,maxlevel=5,side='lower'){
    newlim <- (ll+ul)/2
    fit <- initial2time(x,anames=aname,avalues=newlim,
                        anchor=anchor,type=type,model=model)
    if (fit$LL<(LLmax+LLbuffer)){ # not far enough
        if (maxlevel<1){
            if (side=='lower'){
                return(ll)
            } else {
                return(ul)
            }
        } else {
            if (side=='lower'){
                ul <- newlim
            } else {
                ll <- newlim
            }
        }
    } else { # too far
        if (maxlevel<1){
            return(newlim)
        } else {
            if (side=='lower'){
                ll <- newlim
            } else {
                ul <- newlim
            }
        }        
    }
    recursivelimitsearch_a(aname=aname,ll=ll,ul=ul,LLmax=LLmax,
                           x=x,anchor=anchor,type=type,model=model,
                           maxlevel=maxlevel-1,side=side)
}
getsearchlimits_a <- function(fit,x,anchor=0,type='joint',
                              maxlevel=5,model=1,debug=FALSE){
    if (debug) browser()
    helper <- function(aname,fit,x=x,anchor=0,type=1,model=1,
                       maxlevel=5,alim=c(0,20,1e-5)){
        McL <- mclean(tt=exp(fit$par['t']),d=x$d)
        par <- fit$par
        iname <- paste0(aname,'i')
        par[aname] <- unname(McL[[iname]])
        LLmax <- LL.ludwig(par,x=x,anchor=anchor,type=type,model=model,alim=alim)
        midpoint <- ifelse(McL$truncated,1,unname(McL[[iname]]))
        ll <- recursivelimitsearch_a(aname=aname,ll=alim[1]+alim[3],ul=midpoint,
                                     LLmax=LLmax,x=x,anchor=anchor,type=type,
                                     model=model,maxlevel=maxlevel,side='lower')
        ul <- recursivelimitsearch_a(aname=aname,ll=midpoint,ul=alim[2]-alim[3],
                                     LLmax=LLmax,x=x,anchor=anchor,type=type,
                                     model=model,maxlevel=maxlevel,side='upper')
        c(ll=ll,ul=ul)
    }
    out <- NULL
    if (x$d$U48$option==2 && type%in%c('joint',0,1,3)){
        message('Obtaining U48i search limits')
        lims <- helper(aname='U48',fit=fit,x=x,anchor=anchor,type=type,
                       model=model,maxlevel=maxlevel)
        out <- cbind(out,U48i=lims)
    }
    if (x$d$ThU$option==2 && type%in%c('joint',0,2,4)){
        message('Obtaining ThUi search limits')
        lims <- helper(aname='ThU',fit=fit,x=x,anchor=anchor,
                       type=type,model=model,maxlevel=maxlevel)
        out <- cbind(out,ThUi=lims)
    }
    out
}

time2initial <- function(tt,x=x,type='joint',model=1,alim=c(0,20,1e-5),debug=FALSE){
    if (debug) browser()
    anchor <- c(2,tt)
    init <- init.ludwig(x=x,model=model,anchor=anchor,type=type)
    fit <- stats::optim(init$par,fn=LL.ludwig,method='L-BFGS-B',
                        lower=init$lower,upper=init$upper,hessian=FALSE,
                        x=x,anchor=anchor,type=type,model=model,alim=alim,
                        debug=debug)
    out <- list()
    out$par <- fit$par
    out$LL <- fit$value
    out
}
recursivelimitsearch_t <- function(ll,ul,LLmax,x=x,type=1,maxlevel=5,
                                   model=1,LLbuffer=10,side='lower',
                                   alim=c(0,20,1e-5)){
    if (side=='lower'){
        fit <- time2initial(tt=ll,x=x,type=type,model=model,alim=alim)
    } else {
        fit <- time2initial(tt=ul,x=x,type=type,model=model,alim=alim)
    }
    dl <- (ul-ll)/5
    if (fit$LL<(LLmax+LLbuffer)){ # not far enough
        if (side=='lower'){
            ll <- ll-dl
        } else {
            ul <- ul+dl
        }
    } else {
        if (side=='lower'){
            ll <- ll+dl
        } else {
            ul <- ul-dl
        }
    }
    if (maxlevel<1){
        if (side=='lower'){
            return(ll)
        } else {
            return(ul)
        }
    }
    recursivelimitsearch_t(ll=ll,ul=ul,LLmax=LLmax,
                           x=x,type=type,model=model,
                           maxlevel=maxlevel-1,side=side,alim=alim)
}
getsearchlimits_t <- function(init,x,type='joint',maxlevel=5,
                              model=1,alim=c(0,20,1e-5),debug=FALSE){
    if (debug) browser()
    message('Obtaining t search limits')
    LLmax <- time2initial(tt=init[1],x=x,type=type,model=model,alim=alim)$LL
    ll <- recursivelimitsearch_t(ll=init[2],ul=init[1],LLmax=LLmax,
                                 x=x,type=type,model=model,
                                 maxlevel=maxlevel,side='lower')
    ul <- recursivelimitsearch_t(ll=init[1],ul=init[3],LLmax=LLmax,
                                 x=x,type=type,model=model,
                                 maxlevel=maxlevel,side='upper')
    c(ll,ul)
}

bayeslud <- function(fit,x,anchor=0,type='joint',model=1,
                     alim=c(0,20,1e-5),nsteps=NULL,plot=FALSE){
    if (is.null(nsteps)){
        if (x$d$U48$option==2 && x$d$ThU$option==2) nsteps <- 20
        else nsteps <- 30
    }
    lims <- getsearchlimits_a(fit=fit,x=x,anchor=anchor,type=type,
                              model=model,maxlevel=10)
    pnames <- names(fit$par)
    inames <- colnames(lims)
    np <- length(pnames)
    ilist <- iilist <- list()
    for (iname in inames){
        ilist[[iname]] <- seq(from=lims['ll',iname],
                              to=lims['ul',iname],length.out=nsteps)
        iilist[[iname]] <- 1:nsteps
    }
    igrid <- data.matrix(expand.grid(ilist))
    iigrid <- data.matrix(expand.grid(iilist))
    ng <- nrow(igrid)
    LLgrid <- matrix(NA,ng,np+1)
    colnames(LLgrid) <- c(pnames,'LL')
    aname1 <- ifelse(inames[1]=='U48i','U48','ThU')
    if (length(inames)==1){
        anames <- aname1
    } else {
        aname2 <- ifelse(inames[2]=='U48i','U48','ThU')
        anames <- c(aname1,aname2)
    }
    message('Calculating posterior distribution of the initial activity ratios')
    for (i in 1:ng){
        message('Iteration ',i,'/',ng)
        tfit <- initial2time(x=x,anames=anames,avalues=igrid[i,inames],
                             anchor=anchor,type=type,model=model)
        LLgrid[i,pnames] <- tfit$par[pnames]
        LLgrid[i,inames] <- igrid[i,inames]
        LLgrid[i,'LL'] <- -tfit$LL
    }
    out <- NULL
    for (iname in inames){
        LL <- marginal(LLgrid,iigrid,iilist,iname=iname)
        dx <- diff(ilist[[iname]])
        L <- exp(LL-log_sum_exp(LL+log(c(dx,tail(dx,n=1)))))
        out[[iname]] <- cbind(x=ilist[[iname]],L=L)
    }
    message('Calculating posterior distribution of the age')
    if ('t'%in%pnames){
        init <- exp(c(fit$par['t'],range(LLgrid[,'t'])))
        lims <- getsearchlimits_t(init=init,x=x,type=type,
                                  model=model,maxlevel=10)
        tt <- seq(from=lims[1],to=lims[2],length.out=nsteps)
        init <- lower <- upper <- rep(NA,np-1)
        names(init) <- names(lower) <- names(upper) <- pnames[-1]
        for (pname in names(init)){
            lower[pname] <- min(LLgrid[,pname])
            upper[pname] <- max(LLgrid[,pname])
            if (upper[pname]==lower[pname] &&
                pname%in%c('U48i','ThUi','RaUi','PaUi')){
                lower[pname] <- alim[1]+alim[3]
                upper[pname] <- alim[2]-alim[3]
            }
        }
        LLgridt <- LLgrid[1:nsteps,]
        for (i in 1:nsteps){
            message('Iteration ',i,'/',nsteps)
            for (pname in names(init)){
                init[pname] <- approx(x=exp(LLgrid[,'t']),
                                      y=LLgrid[,pname],xout=tt[i],rule=2)$y
            }
            ifit <- time2initial(tt=tt[i],x=x,type=type,model=model)
            LLgridt[i,'t'] <- log(tt[i])
            LLgridt[i,names(ifit$par)] <- ifit$par
            LLgridt[i,'LL'] <- -ifit$LL
        }
        dt <- diff(tt)
        L <- exp(LLgridt[,'LL']-log_sum_exp(LLgridt[,'LL']+log(c(dt,tail(dt,n=1)))))
        out[['t']] <- cbind(x=tt,L=L)
    }
    if (plot){
        nbpar <- length(out)
        op <- par(mfrow=c(1,nbpar))
        for (bpar in 1:nbpar){
            plot(out[[bpar]],type='b',xlab=names(out)[bpar])
        }
        par(op)
    }
    out
}

marginal <- function(LLgrid,iigrid,iilist,iname='U48i'){
    ii <- iilist[[iname]]
    LL <- rep(NA,length(ii))
    for (i in ii){
        j <- iigrid[,iname]%in%i
        LL[i] <- log_sum_exp(LLgrid[j,'LL'])
    }
    LL
}

prior <- function(x,alim=c(0,20),mean=logit(1,m=alim[1],M=alim[2]),sd=10,log=TRUE){
    lx <- logit(x,m=alim[1],M=alim[2])
    dnorm(lx,mean=mean,sd=sd,log=log)
}

logit <- function(x,m=0,M=1,inverse=FALSE){
    if (inverse){
        out <- m+(M-m)*exp(x)/(1+exp(x))
    } else {
        out <- log(x-m)/log(M-x)
    }
    out
}
