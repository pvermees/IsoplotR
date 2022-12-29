initial2time <- function(x,anames,values,anchor=0,type='joint',model=1,debug=FALSE){
    if (debug){
        browser()
    }
    X <- x
    for (i in seq_along(anames)){
        X$d[[anames[i]]] <- list(x=values[i],sx=0,option=1)
    }
    init <- init.ludwig(x=X,model=model,anchor=anchor,type=type,debug=debug)
    fit <- stats::optim(init$par,fn=LL.ludwig,method='L-BFGS-B',
                       lower=init$lower,upper=init$upper,hessian=FALSE,
                       x=x,X=X,anchor=anchor,type=type,model=model)
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

time2initial <- function(tt,x=x,init,lower,upper,type='joint',model=1){
    anchor <- c(2,tt)
    fit <- stats::optim(init,fn=LL.ludwig,method='L-BFGS-B',
                        lower=lower,upper=upper,hessian=FALSE,
                        x=x,anchor=anchor,type=type,model=model)
    out <- list()
    out$par <- fit$par
    out$LL <- fit$value
    out
}

recursivelimitsearch <- function(aname,ll,ul,LLmax,x=x,anchor=0,
                                 type=1,model=1,buffer=20,maxlevel=5,side='lower'){
    newlim <- (ll+ul)/2
    fit <- initial2time(x,anames=aname,values=newlim,
                        anchor=anchor,type=type,model=model)
    if (fit$LL<(LLmax+buffer)){ # not far enough
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
    recursivelimitsearch(aname=aname,ll=ll,ul=ul,LLmax=LLmax,
                         x=x,anchor=anchor,type=type,model=model,
                         maxlevel=maxlevel-1,side=side)
}

searchlimithelper <- function(aname,fit,x=x,anchor=0,type=1,
                              model=1,maxlevel=5,debug=FALSE){
    if (debug) browser()
    McL <- mclean(tt=exp(fit$par['t']),d=x$d)
    par <- fit$par
    iname <- paste0(aname,'i')
    par[aname] <- unname(McL[[iname]])
    LLmax <- LL.ludwig(par,x=x,anchor=anchor,type=type,model=model)
    midpoint <- ifelse(McL$truncated,1,unname(McL[[iname]]))
    ll <- recursivelimitsearch(aname=aname,ll=0,ul=midpoint,
                               LLmax=LLmax,x=x,anchor=anchor,type=type,
                               model=model,maxlevel=maxlevel,side='lower')
    ul <- recursivelimitsearch(aname=aname,ll=midpoint,ul=20,
                               LLmax=LLmax,x=x,anchor=anchor,type=type,
                               model=model,maxlevel=maxlevel,side='upper')
    c(ll=ll,ul=ul)
}

getsearchlimits <- function(fit,x,anchor=0,type='joint',maxlevel=5,
                            model=1,debug=FALSE){
    out <- NULL
    if (x$d$U48$option==2){
        message('Obtaining U48i search limits')
        lims <- searchlimithelper(aname='U48',fit=fit,x=x,
                                  anchor=anchor,type=type,model=model,
                                  maxlevel=maxlevel,debug=debug)
        out <- cbind(out,U48i=lims)
    }
    if (x$d$ThU$option==2){
        message('Obtaining ThUi search limits')
        lims <- searchlimithelper(aname='ThU',fit=fit,x=x,anchor=anchor,
                                  type=type,model=model,maxlevel=maxlevel)
        out <- cbind(out,ThUi=lims)
    }
    out
}

bayeslud <- function(fit,x,anchor=0,type='joint',model=1,nsteps=NULL,debug=FALSE){
    if (is.null(nsteps)){
        if (x$d$U48$option==2 && x$d$ThU$option==2) nsteps <- 20
        else nsteps <- 30
    }
    lims <- getsearchlimits(fit=fit,x=x,anchor=anchor,type=type,
                            model=model,maxlevel=10,debug=FALSE)
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
        tfit <- initial2time(x=x,anames=anames,values=igrid[i,inames],
                             anchor=anchor,type=type,model=model,debug=FALSE)
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
        mint <- exp(min(LLgrid[,'t']))
        maxt <- exp(max(LLgrid[,'t']))
        tt <- seq(from=mint,to=maxt,length.out=nsteps)
        init <- lower <- upper <- rep(NA,np-1)
        names(init) <- names(lower) <- names(upper) <- pnames[-1]
        for (pname in names(init)){
            lower[pname] <- min(LLgrid[,pname])
            upper[pname] <- max(LLgrid[,pname])
        }
        LLgridt <- LLgrid[1:nsteps,]
        for (i in 1:nsteps){
            message('Iteration ',i,'/',nsteps)
            for (pname in names(init)){
                init[pname] <- approx(x=exp(LLgrid[,'t']),
                                      y=LLgrid[,pname],xout=tt[i])$y
            }
            ifit <- time2initial(tt=tt[i],x=x,init,lower,upper,
                                 type=type,model=model)
            LLgridt[i,'t'] <- log(tt[i])
            LLgridt[i,names(ifit$par)] <- ifit$par
            LLgridt[i,'LL'] <- -ifit$LL
        }
        dt <- diff(tt)
        L <- exp(LLgridt[,'LL']-log_sum_exp(LLgridt[,'LL']+log(c(dt,tail(dt,n=1)))))
        out[['t']] <- cbind(x=tt,L=L)
    }
    if (debug){
        nact <- length(out)
        op <- par(mfrow=c(1,nact))
        for (act in 1:nact){
            plot(out[[act]],type='b')
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
