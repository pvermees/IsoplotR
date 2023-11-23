initial2time <- function(x,anames,avalues,anchor=0,
                         type='joint',model=1){
    X <- x
    for (i in seq_along(anames)){
        X$d[[anames[i]]]$x <- avalues[i]
        X$d[[anames[i]]]$sx <- 0
        X$d[[anames[i]]]$option <- 1
    }
    init <- init.ludwig(x=X,model=model,anchor=anchor,type=type)
    fit <- stats::optim(init$par,fn=LL.ludwig,method='L-BFGS-B',
                        lower=init$lower,upper=init$upper,hessian=FALSE,
                        x=x,X=X,anchor=anchor,type=type,model=model)
    lt <- ifelse('t'%in%names(fit$par),fit$par['t'],anchor[2])
    McL <- mclean(tt=exp(lt),d=X$d)
    out <- list()
    out$par <- fit$par
    for (aname in anames){
        iname <- aname2pname(aname)
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
                           LLbuffer=LLbuffer,x=x,anchor=anchor,type=type,
                           model=model,maxlevel=maxlevel-1,side=side)
}
getsearchlimits_a <- function(fit,x,anchor=0,type='joint',
                              maxlevel=5,LLbuffer=10,model=1){
    helper <- function(aname,fit,x=x,anchor=0,type=1,model=1,maxlevel=5){
        m <- x$d[[aname]]$m
        M <- x$d[[aname]]$M
        buffer <- x$d$buffer
        midpoint <- unname(fit$par[aname2pname(aname)])
        LLmax <- fit$value
        ll <- recursivelimitsearch_a(aname=aname,ll=m+buffer,ul=midpoint,
                                     LLmax=LLmax,LLbuffer=LLbuffer,x=x,
                                     anchor=anchor,type=type,model=model,
                                     maxlevel=maxlevel,side='lower')
        ul <- recursivelimitsearch_a(aname=aname,ll=midpoint,ul=M-buffer,
                                     LLmax=LLmax,LLbuffer=LLbuffer,x=x,
                                     anchor=anchor,type=type,model=model,
                                     maxlevel=maxlevel,side='upper')
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
getsearchlimits_t <- function(LLgrid,LLbuffer,x,fit,type,model){
    message('Obtaining search limits for t')
    mini <- which.min(LLgrid[,'t'])
    maxi <- which.max(LLgrid[,'t'])
    mint <- exp(LLgrid[mini,'t'])
    maxt <- exp(LLgrid[maxi,'t'])
    dt <- maxt-mint
    LLmint <- -LLgrid[mini,'LL']
    LLmax <- fit$value
    lims <- init <- fit$par[-1]
    for (i in 1:10){
        if ((LLmint < (LLmax+LLbuffer)) && (mint > dt/4)){
            mint <- (mint-dt/4)
            ifit <- stats::optim(init,fn=LL.ludwig,hessian=FALSE,x=x,
                                 anchor=c(2,mint),type=type,model=model)
            LLmint <- ifit$value
        } else {
            break
        }
    }
    LLmaxt <- -LLgrid[maxi,'LL']
    for (i in 1:10){
        if (LLmaxt < (LLmax+LLbuffer)){
            maxt <- (maxt+dt/4)
            ifit <- stats::optim(init,fn=LL.ludwig,hessian=FALSE,x=x,
                                 anchor=c(2,maxt),type=type,model=model)
            LLmaxt <- ifit$value
        } else {
            break
        }
    }
    c(mint,maxt)
}

bayeslud <- function(fit,x,anchor=0,type='joint',model=1,
                     nsteps=NULL,plot=FALSE){
    if (is.null(nsteps)){
        if (x$d$U48$option==2 && x$d$ThU$option==2) nsteps <- 20
        else nsteps <- 30
    }
    LLbuffer <- 10
    lims <- getsearchlimits_a(fit=fit,x=x,anchor=anchor,type=type,
                              LLbuffer=LLbuffer,model=model,maxlevel=10)
    inames <- colnames(lims)
    pnames <- names(fit$par)
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
    LLgrid <- matrix(NA,nrow=ng,ncol=np+1)
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
        L <- exp(LL-log_sum_exp(LL+log(c(dx,utils::tail(dx,n=1)))))
        out[[iname]] <- cbind(x=ilist[[iname]],L=L)
    }
    if ('t'%in%pnames){
        lims <- getsearchlimits_t(LLgrid,LLbuffer=LLbuffer,x=x,
                                  fit=fit,type=type,model=model)
        tt <- seq(from=lims[1],to=lims[2],length.out=nsteps)
        ti <- which('t' %in% pnames)
        lower <- upper <- init <- matrix(NA,nrow=nsteps,ncol=np-1)
        colnames(lower) <- colnames(upper) <- colnames(init) <- pnames[-ti]
        for (pname in pnames[-ti]){
            init[,pname] <- stats::approx(x=LLgrid[,'t'],y=LLgrid[,pname],
                                          xout=log(tt),rule=2)$y
        }
        if ('a0' %in% pnames){
            lower[,'a0'] <- init[,'a0']-1
            upper[,'a0'] <- init[,'a0']+1
        }
        if ('b0' %in% pnames){
            lower[,'b0'] <- init[,'b0']-1
            upper[,'b0'] <- init[,'b0']+1
        }
        if ('U48i' %in% pnames){
            lower[,'U48i'] <- init[,'U48i']/2
            upper[,'U48i'] <- init[,'U48i']*2
        }
        if ('ThUi' %in% pnames){
            lower[,'ThUi'] <- init[,'ThUi']/2
            upper[,'ThUi'] <- init[,'ThUi']*2
        }
        LLgridt <- LLgrid[1:nsteps,]
        message('Calculating posterior distribution of the age')
        for (i in 1:nsteps){
            message('Iteration ',i,'/',nsteps)
            anchor <- c(2,tt[i])
            ifit <- contingencyfit(par=init[i,],fn=LL.ludwig,lower=lower[i,],
                                   upper=upper[i,],hessian=FALSE,x=x,
                                   anchor=anchor,type=type,model=model)
            LLgridt[i,'t'] <- tt[i]
            LLgridt[i,names(ifit$par)] <- ifit$par
            LLgridt[i,'LL'] <- -ifit$value
        }
        L <- exp(LLgridt[,'LL'] - log_sum_exp(LLgridt[,'LL']))
        out[['t']] <- cbind(x=tt,L=L)
    }
    if (plot){
        nbpar <- length(out)
        op <- graphics::par(mfrow=c(1,nbpar))
        for (bpar in names(out)){
            plot(out[[bpar]],type='b',xlab=bpar)
            if (bpar=='t') xx <- exp(fit$par[bpar])
            else xx <- fit$par[bpar]
            graphics::lines(rep(xx,2),range(out[[bpar]][,2]))
        }
        graphics::par(op)
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

prior <- function(x,a,log=TRUE){
    lx <- logit(x,m=a$m,M=a$M)
    mu <- logit(a$x0,m=a$m,M=a$M)
    stats::dnorm(lx,mean=mu,sd=a$sd,log=log)
}

logit <- function(x,m=0,M=1,inverse=FALSE){
    if (inverse){
        out <- m+(M-m)*exp(x)/(1+exp(x))
    } else {
        out <- log(x-m)/log(M-x)
    }
    out
}
