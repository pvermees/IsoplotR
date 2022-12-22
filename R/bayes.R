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
    tt <- ifelse('t'%in%names(fit$par),fit$par['t'],anchor[2])
    McL <- mclean(tt=tt,d=X$d)
    out <- list()
    out$par <- fit$par
    for (aname in anames){
        iname <- paste0(aname,'i')
        out$par[iname] <- log(McL[[iname]])
    }
    out$LL <- LL.ludwig(out$par,x=x,anchor=anchor,type=type,model=model)
    out
}

recursivelimitsearch <- function(aname,iname,ll,ul,LLmax,method,
                                 x=x,anchor=0,type=1,model=1,
                                 maxlevel=5,side='lower',debug=FALSE){
    newlim <- (ll+ul)/2
    fit <- initial2time(x,anames=aname,values=newlim,method=method,
                        anchor=anchor,type=type,model=model)
    if (fit$LL<(LLmax+100)){ # not far enough
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
    recursivelimitsearch(aname=aname,iname=iname,ll=ll,ul=ul,LLmax=LLmax,
                         method=method,x=x,anchor=anchor,type=type,
                         model=model,maxlevel=maxlevel-1,side=side)
}

searchlimithelper <- function(aname,iname,fit,method,x=x,
                              anchor=0,type=1,model=1,debug=FALSE){
    McL <- mclean(tt=exp(fit$par['t']),d=x$d)
    par <- fit$par
    par[aname] <- log(unname(McL[[iname]]))
    LLmax <- LL.ludwig(par,x=x,anchor=anchor,type=type,model=model)
    ll <- recursivelimitsearch(aname=aname,iname=iname,ll=0,ul=unname(McL[[aname]]),
                               LLmax=LLmax,method=method,x=x,anchor=anchor,
                               type=type,model=model,maxlevel=5,side='lower')
    ul <- recursivelimitsearch(aname=aname,iname=iname,ll=unname(McL[[aname]]),ul=20,
                               LLmax=LLmax,method=method,x=x,anchor=anchor,
                               type=type,model=model,maxlevel=5,side='upper')
    c(ll=ll,ul=ul)
}

getsearchlimits <- function(fit,x,anchor=0,method='Nelder=-Mead',
                            type='joint',model=1){
    out <- NULL
    if (x$d$U48$option==2){
        lims <- searchlimithelper(lims,aname='U48',iname='U48i',fit=fit,
                                  method=method,x=x,anchor=anchor,
                                  type=type,model=model)
        out <- cbind(out,U48i=lims)
    }
    if (x$d$ThU$option==2){
        lims <- searchlimithelper(lims,aname='ThU',iname='ThUi',fit=fit,
                                  method=method,x=x,anchor=anchor,
                                  type=type,model=model)
        out <- cbind(out,ThUi=lims)
    }
    out
}

bayeslud <- function(fit,x,anchor=0,type='joint',model=1,debug=FALSE){
    if (length(fit$par)>1) method <- "Nelder-Mead"
    else method <- "BFGS"
    lims <- getsearchlimits(fit=fit,x=x,method=method,
                            anchor=anchor,type=type,model=model)
    pnames <- names(fit$par)
    inames <- colnames(lims)
    ni <- length(inames)
    np <- length(pnames)
    nn <- 20
    ilist <- list()
    for (iname in inames){
        ilist[[iname]] <- seq(from=lims['ll',iname],
                                to=lims['ul',iname],
                                length.out=nn)
    }
    pargrid <- expand.grid(ilist)
    ng <- nrow(pargrid)
    out <- matrix(NA,ng,np+ni+1)
    colnames(out) <- c(pnames,inames,'LL')
    aname1 <- ifelse(iname[1]=='U48i','U48','ThU')
    if (ni==1){
        anames <- aname1
    } else {
        aname2 <- ifelse(iname[2]=='U48i','U48','ThU')
        anames <- c(aname1,aname2)
    }
    for (i in 1:ng){
        fit <- initial2time(x=x,anames=anames,values=pargrid[i,inames],
                            method=method,anchor=anchor,type=type,
                            model=model,debug=FALSE)
        out[i,pnames] <- fit$par[pnames]
        out[i,inames] <- pargrid[i,inames]
        out[i,'LL'] <- fit$LL
    }
    if (debug){
        L <- exp(min(out[,'LL'])-out[,'LL'])
        plot(out[,'U48i'],L,type='b')
    }
    out
}
