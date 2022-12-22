initial2time <- function(x,anames,values,method='Nelder-Mead',
                         anchor=0,type=1,model=1,debug=FALSE){
    if (debug){
        browser()
    }
    X <- x
    for (i in seq_along(anames)){
        X$d[[anames[i]]] <- list(x=values[i],sx=0,option=1)
    }
    init <- init.ludwig(x=X,model=model,anchor=anchor,type=type)
    fit <- stats::optim(init,fn=LL.ludwig,method=method,hessian=FALSE,
                        x=X,anchor=anchor,type=type,model=model)
    lt <- ifelse('t'%in%names(fit$par),fit$par['t'],anchor[2])
    McL <- mclean(tt=exp(lt),d=X$d)
    out <- list()
    out$par <- fit$par
    for (aname in anames){
        iname <- paste0(aname,'i')
        out$par[iname] <- log(McL[[iname]])
    }
    out$LL <- LL.ludwig(out$par,x=x,anchor=anchor,type=type,model=model)
    out
}

recursivelimitsearch <- function(aname,ll,ul,LLmax,method,
                                 x=x,anchor=0,type=1,model=1,
                                 maxlevel=5,side='lower'){
    newlim <- (ll+ul)/2
    fit <- initial2time(x,anames=aname,values=newlim,method=method,
                        anchor=anchor,type=type,model=model)
    if (fit$LL<(LLmax+25)){ # not far enough
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
                         method=method,x=x,anchor=anchor,type=type,
                         model=model,maxlevel=maxlevel-1,side=side)
}

searchlimithelper <- function(aname,fit,method,x=x,anchor=0,
                              type=1,model=1,maxlevel=5){
    McL <- mclean(tt=exp(fit$par['t']),d=x$d)
    par <- fit$par
    iname <- paste0(aname,'i')
    par[aname] <- log(unname(McL[[iname]]))
    LLmax <- LL.ludwig(par,x=x,anchor=anchor,type=type,model=model)
    ll <- recursivelimitsearch(aname=aname,ll=0,ul=unname(McL[[iname]]),
                               LLmax=LLmax,method=method,x=x,anchor=anchor,
                               type=type,model=model,maxlevel=maxlevel,side='lower')
    ul <- recursivelimitsearch(aname=aname,ll=unname(McL[[iname]]),ul=20,
                               LLmax=LLmax,method=method,x=x,anchor=anchor,
                               type=type,model=model,maxlevel=maxlevel,side='upper')
    c(ll=ll,ul=ul)
}

getsearchlimits <- function(fit,x,anchor=0,method='Nelder=-Mead',
                            type='joint',maxlevel=5,model=1){
    out <- NULL
    if (x$d$U48$option==2){
        message('Obtaining U48i search limits')
        lims <- searchlimithelper(aname='U48',fit=fit,
                                  method=method,x=x,anchor=anchor,
                                  type=type,model=model,maxlevel=maxlevel)
        out <- cbind(out,U48i=lims)
    }
    if (x$d$ThU$option==2){
        message('Obtaining ThUi search limits')
        lims <- searchlimithelper(aname='ThU',fit=fit,
                                  method=method,x=x,anchor=anchor,
                                  type=type,model=model,maxlevel=maxlevel)
        out <- cbind(out,ThUi=lims)
    }
    out
}

bayeslud <- function(fit,x,anchor=0,type='joint',model=1,debug=FALSE,nsteps=15){
    if (length(fit$par)>1) method <- "Nelder-Mead"
    else method <- "BFGS"
    lims <- getsearchlimits(fit=fit,x=x,method=method,anchor=anchor,
                            type=type,model=model,maxlevel=10)
    pnames <- names(fit$par)
    inames <- colnames(lims)
    ni <- length(inames)
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
    LLgrid <- matrix(NA,ng,np+ni+1)
    colnames(LLgrid) <- c(pnames,inames,'LL')
    aname1 <- ifelse(inames[1]=='U48i','U48','ThU')
    if (ni==1){
        anames <- aname1
    } else {
        aname2 <- ifelse(inames[2]=='U48i','U48','ThU')
        anames <- c(aname1,aname2)
    }
    for (i in 1:ng){
        message('Iteration ',i,'/',ng)
        fit <- initial2time(x=x,anames=anames,
                            values=igrid[i,inames],
                            method=method,anchor=anchor,type=type,
                            model=model,debug=FALSE)
        LLgrid[i,pnames] <- fit$par[pnames]
        LLgrid[i,inames] <- igrid[i,inames]
        LLgrid[i,'LL'] <- -fit$LL
    }
    out <- NULL
    for (iname in inames){
        x <- ilist[[iname]]
        LL <- marginal(LLgrid,iigrid,iilist,iname=iname)
        out[[iname]] <- cbind(x=x,LL=LL)
    }
    if (debug){
        xLL <- out[['U48i']]
        dx <- diff(xLL[,'x'])
        y <- exp(xLL[,'LL'])/sum(exp(xLL[,'LL']*c(dx,tail(dx,n=1))))
        plot(x,y,type='b')
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
