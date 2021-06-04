# returns the lower and upper intercept age (for Wetherill concordia)
# or the lower intercept age and 207Pb/206Pb intercept (for Tera-Wasserburg)
concordia.intersection.ludwig <- function(x,wetherill=TRUE,exterr=FALSE,
                                          alpha=0.05,model=1,anchor=0){
    fit <- ludwig(x,exterr=exterr,model=model,anchor=anchor)
    out <- fit
    out$fact <- tfact(alpha,fit$df)
    out$format <- x$format
    if (wetherill & !measured.disequilibrium(x$d)){
        wfit <- twfit2wfit(fit,x)
        out$par <- wfit$par
        out$cov <- wfit$cov
        names(out$par) <- c('t[l]','t[u]')
    } else {
        out$par <- fit$par
        out$cov <- fit$cov
    }
    np <- length(out$par)
    out$alpha <- alpha
    if (inflate(out)){
        out$err <- matrix(NA,3,np)
        rownames(out$err) <- c('s','ci','disp')
        out$err['disp',] <-
            out$fact*sqrt(fit$mswd)*sqrt(diag(out$cov))
    } else {
        out$err <- matrix(NA,2,np)
        rownames(out$err) <- c('s','ci')
    }
    if (model==3) out$fact <- nfact(alpha)
    out$err['s',] <- sqrt(diag(out$cov))
    out$err['ci',] <- out$fact*out$err['s',]
    colnames(out$err) <- names(out$par)
    out
}
# extracts concordia intersection parameters from an ordinary York fit
concordia.intersection.ab <- function(a,b,covmat=matrix(0,2,2),
                                      exterr=FALSE,wetherill=FALSE,d=diseq()){
    l8 <- lambda('U238')[1]
    ta <- get.Pb207Pb206.age(a,d=d)[1]
    out <- c(1,a) # tl, 7/6 intercept
    if (wetherill) names(out) <- c('t[l]','t[u]')
    else names(out) <- c('t[l]','76i')
    if (b<0) { # negative slope => two (or zero) intersections with concordia line
        tb <- get.Pb206U238.age(-b/a,d=d)[1]
        tlu <- recursive.search(tm=tb,tM=ta,a=a,b=b,d=d)
        out['t[l]'] <- tlu[1]
        if (wetherill) out['t[u]'] <- tlu[2]
    } else {
        search.range <- c(ta,2/l8)
        out['t[l]'] <- stats::uniroot(intersection.misfit.york,
                                      interval=search.range,a=a,b=b,d=d)$root
    }
    out
}

recursive.search <- function(tm,tM,a,b,d=diseq(),depth=1){
    out <- c(NA,NA)
    if (depth<3){
        mid <- (tm+tM)/2
        mfmin <- intersection.misfit.york(tm,a=a,b=b,d=d)
        mfmid <- intersection.misfit.york(mid,a=a,b=b,d=d)
        mfmax <- intersection.misfit.york(tM,a=a,b=b,d=d)
        if (mfmin*mfmid<0){ # different signs
            out[1] <- stats::uniroot(intersection.misfit.york,
                                     interval=c(tm,mid),a=a,b=b,d=d)$root
        } else {
            out <- recursive.search(tm=tm,tM=mid,a=a,b=b,d=d,depth=depth+1)
        }
        if (mfmax*mfmid<0){ # different signs
            out[2] <- stats::uniroot(intersection.misfit.york,
                                     interval=c(mid,tM),a=a,b=b,d=d)$root
        } else {
            tlu <- recursive.search(tm=mid,tM=tM,a=a,b=b,d=d,depth=depth+1)
            if (is.na(out[1])) out[1] <- tlu[1]
            if (is.na(out[2])) out[2] <- tlu[2]
        }
        if (all(is.na(out))){ # no intersection
            tlu <- stats::optimise(intersection.misfit.york,
                                   interval=c(tm,tM),a=a,b=b,d=d)$minimum
            out <- rep(tlu,2)
        }
    }
    out
}

# extract the lower and upper discordia intercept from the parameters
# of a Ludwig fit (initial Pb ratio and lower intercept age)
twfit2wfit <- function(fit,x){
    tt <- fit$par['t']
    buffer <- 1 # start searching 1Ma above or below first intercept age
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    U <- iratio('U238U235')[1]
    E <- matrix(0,3,3)
    J <- matrix(0,2,3)
    if (x$format %in% c(1,2,3)){
        a0 <- 1
        b0 <- fit$par['76i']
        E[c(1,3),c(1,3)] <- fit$cov[1:2,1:2]
    } else if (x$format %in% c(4,5,6)){
        a0 <- fit$par['64i']
        b0 <- fit$par['74i']
        E <- fit$cov[1:3,1:3]
    } else if (x$format%in%c(7,8)){
        a0 <- fit$par['68i']
        b0 <- fit$par['78i']
        E <- fit$cov[1:3,1:3]
    } else {
        stop('Incorrect input format')
    }
    disc.slope <- a0/(b0*U)
    conc.slope <- (l8*exp(l8*tt))/(l5*exp(l5*tt))
    if (disc.slope < conc.slope){
        search.range <- c(tt,get.Pb207Pb206.age(b0/a0,d=x$d)[1])+buffer
        tl <- tt
        tu <- stats::uniroot(intersection.misfit.ludwig,interval=search.range,
                             t2=tt,a0=a0,b0=b0,d=x$d)$root
    } else if (disc.slope < l8/l5){
        search.range <- c(-1000,tt-buffer)
        tl <- stats::uniroot(intersection.misfit.ludwig,interval=search.range,
                             t2=tt,a0=a0,b0=b0,d=x$d)$root
        tu <- tt
    } else { # only one intercept
        tl <- -1000
        tu <- tt
    }
    du <- mclean(tt=tu,d=x$d)
    dl <- mclean(tt=tl,d=x$d)
    XX <- du$Pb207U235 - dl$Pb207U235
    YY <- du$Pb206U238 - dl$Pb206U238
    BB <- a0/(b0*U)
    D <- (YY-BB*XX)^2 # misfit
    dXX.dtu <-  du$dPb207U235dt
    dXX.dtl <- -dl$dPb207U235dt
    dYY.dtu <-  du$dPb206U238dt
    dYY.dtl <- -dl$dPb206U238dt
    dBB.da0 <-  1/(b0*U)
    dBB.db0 <- -BB/b0
    dD.dtl <- 2*(YY-BB*XX)*(dYY.dtl-BB*dXX.dtl)
    dD.dtu <- 2*(YY-BB*XX)*(dYY.dtu-BB*dXX.dtu)
    dD.da0 <- 2*(YY-BB*XX)*(-dBB.da0*XX)
    dD.db0 <- 2*(YY-BB*XX)*(-dBB.db0*XX)
    if (conc.slope > disc.slope){
        J[1,1] <- 1
        J[2,1] <- -dD.dtl/dD.dtu
        J[2,2] <- -dD.da0/dD.dtu
        J[2,3] <- -dD.db0/dD.dtu
    } else {
        J[1,1] <- -dD.dtu/dD.dtl
        J[1,2] <- -dD.da0/dD.dtl
        J[1,3] <- -dD.db0/dD.dtl
        J[2,1] <- 1
    }
    out <- list()
    out$par <- c(tl,tu)
    out$cov <- J %*% E %*% t(J)
    out
}

# t1 = 1st Wetherill intercept, t2 = 2nd Wetherill intercept
# a0 = 64i, b0 = 74i on TW concordia
intersection.misfit.ludwig <- function(t1,t2,a0,b0,d=diseq()){
    tl <- min(t1,t2)
    tu <- max(t1,t2)
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    U <- iratio('U238U235')[1]
    du <- mclean(tt=tu,d=d)
    dl <- mclean(tt=tl,d=d)
    XX <- du$Pb207U235 - dl$Pb207U235
    YY <- du$Pb206U238 - dl$Pb206U238
    BB <- a0/(b0*U)
    # misfit is based on difference in slope in Wetherill space
    YY - BB*XX
}
# a = intercept, b = slope on TW concordia
intersection.misfit.york <- function(tt,a,b,d=diseq()){
    D <- mclean(tt=tt,d=d)
    # misfit is based on difference in slope in TW space
    #D$Pb207U235/U - a*D$Pb206U238 - b
    (D$Pb207Pb206-a)*D$Pb206U238 - b
}

discordia.line <- function(fit,wetherill,d=diseq()){
    X <- c(0,0)
    Y <- c(0,0)
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    J <- matrix(0,1,2)
    usr <- graphics::par('usr')
    if (wetherill){
        if (measured.disequilibrium(d)){
            U85 <- iratio('U238U235')[1]
            fit2d <- tw3d2d(fit)
            xy1 <- age_to_wetherill_ratios(fit$par[1],d=d)
            x1 <- xy1$x[1]
            x2 <- usr[2]
            y1 <- xy1$x[2]
            dydx <- 1/(U85*fit$par[2])
            y2 <- y1 + (x2-x1)*dydx
            X <- c(x1,x2)
            Y <- c(y1,y2)
            cix <- NA # computing confidence envelopes is very tricky
            ciy <- NA # for this rarely used function -> don't bother
        } else {
            tl <- fit$par[1]
            tu <- fit$par[2]
            X <- age_to_Pb207U235_ratio(fit$par,d=d)[,'75']
            Y <- age_to_Pb206U238_ratio(fit$par,d=d)[,'68']
            x <- seq(from=max(0,usr[1],X[1]),to=min(usr[2],X[2]),length.out=50)
            du <- mclean(tt=tu,d=d)
            dl <- mclean(tt=tl,d=d)
            aa <- du$Pb206U238 - dl$Pb206U238
            bb <- x - dl$Pb207U235
            cc <- du$Pb207U235 - dl$Pb207U235
            dd <- dl$Pb206U238
            y <- aa*bb/cc + dd
            dadtl <- -dl$dPb206U238dt
            dbdtl <- -dl$dPb207U235dt
            dcdtl <- -dl$dPb207U235dt
            dddtl <- dl$dPb206U238dt
            dadtu <- du$dPb206U238dt
            dbdtu <- 0
            dcdtu <- du$dPb207U235dt
            dddtu <- 0
            J1 <- dadtl*bb/cc + dbdtl*aa/cc - dcdtl*aa*bb/cc^2 + dddtl # dydtl
            J2 <- dadtu*bb/cc + dbdtu*aa/cc - dcdtu*aa*bb/cc^2 + dddtu # dydtu
            E11 <- fit$cov[1,1]
            E12 <- fit$cov[1,2]
            E22 <- fit$cov[2,2]
            sy <- errorprop1x2(J1,J2,fit$cov[1,1],fit$cov[2,2],fit$cov[1,2])
            ul <- y + fit$fact*sy
            ll <- y - fit$fact*sy
            t75 <- get.Pb207U235.age(x,d=d)[,'t75']
            yconc <- age_to_Pb206U238_ratio(t75,d=d)[,'68']
            overshot <- ul>yconc
            ul[overshot] <- yconc[overshot]
            cix <- c(x,rev(x))
            ciy <- c(ll,rev(ul))
        }
    } else {
        fit2d <- tw3d2d(fit)
        X[1] <- age_to_U238Pb206_ratio(fit2d$par['t'],d=d)[,'86']
        Y[1] <- age_to_Pb207Pb206_ratio(fit2d$par['t'],d=d)[,'76']
        r75 <- age_to_Pb207U235_ratio(fit2d$par['t'],d=d)[,'75']
        r68 <- 1/X[1]
        Y[2] <- fit2d$par['76i']
        xl <- X[1]
        yl <- Y[1]
        y0 <- Y[2]
        tl <- check.zero.UPb(fit2d$par['t'])
        U <- settings('iratio','U238U235')[1]
        nsteps <- 100
        x <- seq(from=max(.Machine$double.xmin,usr[1]),to=usr[2],length.out=nsteps)
        y <- yl + (y0-yl)*(1-x*r68) # = y0 + yl*x*r68 - y0*x*r68
        D <- mclean(tt=tl,d=d)
        d75dtl <- D$dPb207U235dt
        d68dtl <- D$dPb206U238dt
        dyldtl <- (d75dtl*r68 - r75*d68dtl)/(U*r68^2)
        J1 <- dyldtl*x*r68 + yl*x*d68dtl - y0*x*d68dtl # dy/dtl
        J2 <- 1 - x*r68                                # dy/dy0
        sy <- errorprop1x2(J1,J2,fit2d$cov[1,1],fit2d$cov[2,2],fit2d$cov[1,2])
        ul <- y + fit2d$fact*sy
        ll <- y - fit2d$fact*sy
        yconc <- rep(0,nsteps)
        t68 <- get.Pb206U238.age(1/x,d=d)[,'t68']
        yconc <- age_to_Pb207Pb206_ratio(t68,d=d)[,'76']
        # correct overshot confidence intervals:
        if (y0>yl){ # negative slope
            overshot <- (ll<yconc & ll<y0/2)
            ll[overshot] <- yconc[overshot]
            overshot <- (ul<yconc & ul<y0/2)
            ul[overshot] <- yconc[overshot]
        } else {    # positive slope
            overshot <- ul>yconc
            ul[overshot] <- yconc[overshot]
            overshot <- ll>yconc
            ll[overshot] <- yconc[overshot]
        }
        cix <- c(x,rev(x))
        ciy <- c(ll,rev(ul))
    }
    graphics::polygon(cix,ciy,col='gray80',border=NA)
    graphics::lines(X,Y)
}

tw3d2d <- function(fit){
    out <- list(par=fit$par,cov=fit$cov,fact=fit$fact)
    if (fit$format > 3){
        labels <- c('t','76i')
        out$par <- c(fit$par['t'],fit$par[3]/fit$par[2]) # par = c(Pb206i,Pb207i)
        J <- matrix(0,2,3)
        J[1,1] <- 1
        J[2,2] <- -fit$par[3]/fit$par[2]^2
        J[2,3] <- 1/fit$par[2]
        out$cov <- J %*% fit$cov[1:3,1:3] %*% t(J)
        names(out$par) <- labels
        colnames(out$cov) <- labels
    }
    out
}

# this would be much easier in unicode but that doesn't render in PDF:
discordia.title <- function(fit,wetherill,sigdig=2,...){
    lower.age <- roundit(fit$par[1],fit$err[,1],sigdig=sigdig,text=TRUE)
    if (inflate(fit)){
        args1 <- quote(a%+-%b~'|'~c~'|'~d~u~'(n='*n*')')
        args2 <- quote(a%+-%b~'|'~c~'|'~d~u)
    } else {
        args1 <- quote(a%+-%b~'|'~c~u~'(n='*n*')')
        args2 <- quote(a%+-%b~'|'~c~u)
    }
    list1 <- list(a=lower.age[1],b=lower.age[2],
                  c=lower.age[3],u='Ma',n=fit$n)
    if (wetherill){
        upper.age <- roundit(fit$par[2],fit$err[,2],sigdig=sigdig,text=TRUE)
        expr1 <- quote('lower intercept =')
        expr2 <- quote('upper intercept =')
        list2 <- list(a=upper.age[1],b=upper.age[2],c=upper.age[3],u='Ma')
        if (inflate(fit)){
            list1$d <- lower.age[4]
            list2$d <- upper.age[4]
        }
    } else if (fit$format%in%c(1,2,3)){
        i76 <- roundit(fit$par['76i'],fit$err[,'76i'],sigdig=sigdig,text=TRUE)
        expr1 <- quote('age =')
        expr2 <- quote('('^207*'Pb/'^206*'Pb)'[o]*'=')
        list2 <- list(a=i76[1],b=i76[2],c=i76[3],u='')
        if (inflate(fit)){
            list1$d <- lower.age[4]
            list2$d <- i76[4]
        }
    } else if (fit$format%in%c(4,5,6)){
        i64 <- roundit(fit$par['64i'],fit$err[,'64i'],sigdig=sigdig,text=TRUE)
        i74 <- roundit(fit$par['74i'],fit$err[,'74i'],sigdig=sigdig,text=TRUE)
        expr1 <- quote('age =')
        expr2 <- quote('('^206*'Pb/'^204*'Pb)'[o]*'=')
        expr3 <- quote('('^207*'Pb/'^204*'Pb)'[o]*'=')
        list2 <- list(a=i64[1],b=i64[2],c=i64[3],u='')
        list3 <- list(a=i74[1],b=i74[2],c=i74[3],u='')
        if (inflate(fit)){
            list1$d <- lower.age[4]
            list2$d <- i64[4]
            list3$d <- i74[4]
        }
        call3 <- substitute(e~a,list(e=expr3,a=args2))
        line3 <- do.call('substitute',list(call3,list3))        
    } else if (fit$format%in%c(7,8)){
        i86 <- 1/fit$par['68i']
        i87 <- 1/fit$par['78i']
        i86err <- i86*fit$err[,'68i']/fit$par['68i']
        i87err <- i87*fit$err[,'78i']/fit$par['78i']
        ri86 <- roundit(i86,i86err,sigdig=sigdig,text=TRUE)
        ri87 <- roundit(i87,i87err,sigdig=sigdig,text=TRUE)
        expr1 <- quote('age =')
        expr2 <- quote('('^208*'Pb/'^206*'Pb)'[o]*'=')
        expr3 <- quote('('^208*'Pb/'^207*'Pb)'[o]*'=')
        list2 <- list(a=ri86[1],b=ri86[2],c=ri86[3],u='')
        list3 <- list(a=ri87[1],b=ri87[2],c=ri87[3],u='')
        if (inflate(fit)){
            list1$d <- lower.age[4]
            list2$d <- ri86[4]
            list3$d <- ri87[4]
        }
        call3 <- substitute(e~a,list(e=expr3,a=args2))
        line3 <- do.call('substitute',list(call3,list3))
    }
    call1 <- substitute(e~a,list(e=expr1,a=args1))
    call2 <- substitute(e~a,list(e=expr2,a=args2))
    line1 <- do.call('substitute',list(call1,list1))
    line2 <- do.call('substitute',list(call2,list2))
    if (fit$model==1){
        line4 <- substitute('MSWD ='~a*', p('*chi^2*') ='~b,
                            list(a=signif(fit$mswd,sigdig),
                                 b=signif(fit$p.value,sigdig)))
    } else if (fit$model==3){
        ci <- ci_log2lin_lud(fit=fit,fact=fit$fact)
        rounded.disp <- roundit(ci[1],ci[2:3],sigdig=sigdig,text=TRUE)
        line4 <- substitute('overdispersion ='~a+b/-c~'Ma',
                            list(a=rounded.disp[1],b=rounded.disp[3],
                                 c=rounded.disp[2]))
    }
    extrarow <- fit$format>3 & !wetherill
    if (fit$model==1 & extrarow){
        mymtext(line1,line=3,...)
        mymtext(line2,line=2,...)
        mymtext(line3,line=1,...)
        mymtext(line4,line=0,...)
    } else if (fit$model==2 & extrarow){
        mymtext(line1,line=2,...)
        mymtext(line2,line=1,...)
        mymtext(line3,line=0,...)
    } else if (fit$model==3 & extrarow){
        mymtext(line1,line=3,...)
        mymtext(line2,line=2,...)
        mymtext(line3,line=1,...)
        mymtext(line4,line=0,...)
    } else if (fit$model==1){
        mymtext(line1,line=2,...)
        mymtext(line2,line=1,...)
        mymtext(line4,line=0,...)
    } else if (fit$model==2){
        mymtext(line1,line=1,...)
        mymtext(line2,line=0,...)
    } else if (fit$model==3){
        mymtext(line1,line=2,...)
        mymtext(line2,line=1,...)
        mymtext(line4,line=0,...)
    }
}
