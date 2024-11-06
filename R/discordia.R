# returns the lower and upper intercept age (for Wetherill concordia)
# or the lower intercept age and 207Pb/206Pb intercept (for Tera-Wasserburg)
discordia <- function(x,fit,wetherill=TRUE){
    out <- fit
    out$format <- x$format
    if (wetherill){
        wfit <- twfit2wfit(fit,x)
        out$par <- wfit$par
        out$cov <- wfit$cov
    } else {
        out$par <- fit$par
        out$cov <- fit$cov
    }
    np <- length(out$par)
    if (inflate(out)){
        out$err <- matrix(NA,2,np)
        rownames(out$err) <- c('s','disp')
        out$err['disp',] <- sqrt(fit$mswd*diag(out$cov))
    } else {
        out$err <- matrix(NA,1,np)
        rownames(out$err) <- 's'
    }
    colnames(out$err) <- names(out$par)
    out$err['s',] <- sqrt(diag(out$cov))
    out
}

# extract the lower and upper discordia intercept from the parameters
# of a Ludwig fit (initial Pb ratio and lower intercept age)
twfit2wfit <- function(fit,x){
    if (measured_disequilibrium(x$d)) x <- measured2initial(x,fit)
    tt <- fit$par['t']
    buffer <- 1 # start searching 1Ma above or below first intercept age
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    U <- iratio('U238U235')[1]
    if (x$format<4){
        Pb76 <- fit$par['a0']
        dPb76da0 <- 1
        dPb76db0 <- 0
    } else {
        Pb76 <- fit$par['b0']/fit$par['a0']
        dPb76da0 <- -Pb76/fit$par['a0']
        dPb76db0 <- 1/fit$par['a0']
    }
    E <- fit$cov[c('t','a0','b0','w'),c('t','a0','b0','w')]
    J <- matrix(0,3,4)
    md <- mediand(x$d)
    D <- mclean(tt,d=md)
    disc.slope <- 1/(U*Pb76)
    conc.slope <- D$dPb206U238dt/D$dPb207U235dt
    if (disc.slope < conc.slope){
        search.range <- c(tt,get_Pb207Pb206_age(Pb76,d=md)[1])+buffer
        tl <- tt
        tu <- stats::uniroot(intersection_misfit_ludwig,interval=search.range,
                             t2=tt,disc.slope=disc.slope,d=md)$root
    } else {
        search.range <- c(0,tt-buffer)
        if (check_equilibrium(d=x$d)) search.range[1] <- -1000
        tl <- tryCatch(
            stats::uniroot(intersection_misfit_ludwig,
                           interval=search.range,
                           t2=tt,disc.slope=disc.slope,d=md)$root
          , error=function(e){
              stop("Can't find the lower intercept.",
                   "Try fitting the data in Tera-Wasserburg space.")
          })
        tu <- tt
    }
    du <- mclean(tt=tu,d=md)
    dl <- mclean(tt=tl,d=md)
    XX <- du$Pb207U235 - dl$Pb207U235
    YY <- du$Pb206U238 - dl$Pb206U238
    BB <- disc.slope
    D <- (YY-BB*XX)^2 # misfit
    dXX.dtu <-  du$dPb207U235dt
    dXX.dtl <- -dl$dPb207U235dt
    dYY.dtu <-  du$dPb206U238dt
    dYY.dtl <- -dl$dPb206U238dt
    dBB.da0 <-  -dPb76da0*disc.slope/Pb76
    dBB.db0 <- -dPb76db0*disc.slope/Pb76
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
    J[3,4] <- 1
    out <- list()
    out$par <- c(tl,tu,fit$par[-1])
    tokeep <- names(fit$par)[-1]
    pnames <- c('t[l]','t[u]',tokeep)
    np <- length(pnames)
    out$cov <- matrix(0,np,np)
    rownames(out$cov) <- colnames(out$cov) <- names(out$par) <- pnames
    out$cov[tokeep,tokeep] <- fit$cov[tokeep,tokeep]
    toreplace <- c('t[l]','t[u]','w')
    out$cov[toreplace,toreplace] <- J %*% E %*% t(J)
    out
}

# t1 = 1st Wetherill intercept, t2 = 2nd Wetherill intercept
intersection_misfit_ludwig <- function(t1,t2,disc.slope,d=diseq()){
    tl <- min(t1,t2)
    tu <- max(t1,t2)
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    U <- iratio('U238U235')[1]
    du <- mclean(tt=tu,d=d)
    dl <- mclean(tt=tl,d=d)
    XX <- du$Pb207U235 - dl$Pb207U235
    YY <- du$Pb206U238 - dl$Pb206U238
    # misfit is based on difference in slope in Wetherill space
    YY - disc.slope*XX
}
# a = intercept, b = slope on TW concordia
intersection_misfit_york <- function(tt,a,b,d=diseq()){
    D <- mclean(tt=tt,d=d)
    # misfit is based on difference in slope in TW space
    (D$Pb207Pb206-a)*D$Pb206U238 - b
}

discordia_line <- function(fit,wetherill,d=diseq(),oerr=3){
    X <- c(0,0)
    Y <- c(0,0)
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    J <- matrix(0,1,2)
    usr <- graphics::par('usr')
    if (wetherill){
        if (measured_disequilibrium(d)){
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
            X <- age_to_Pb207U235_ratio(c(tl,tu),d=d)[,'75']
            Y <- age_to_Pb206U238_ratio(c(tl,tu),d=d)[,'68']
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
            vy <- errorprop1x2(J1,J2,fit$cov[1,1],fit$cov[2,2],fit$cov[1,2])
            ciy <- ci(x=y,sx=sqrt(vy),oerr=oerr,absolute=TRUE)
            ul <- y + ciy
            ll <- y - ciy
            t75 <- get_Pb207U235_age(x,d=d)[,'t75']
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
        Y[2] <- fit2d$par['a0']
        xl <- X[1]
        yl <- Y[1]
        y0 <- Y[2]
        tl <- check_zero_UPb(fit2d$par['t'])
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
        vy <- errorprop1x2(J1,J2,fit2d$cov[1,1],fit2d$cov[2,2],fit2d$cov[1,2])
        ciy <- ci(x=y,sx=sqrt(vy),oerr=oerr,absolute=TRUE)
        ul <- y + ciy
        ll <- y - ciy
        t68 <- get_Pb206U238_age(1/x,d=d)[,'t68']
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
    out <- list(par=fit$par,cov=fit$cov)
    if (fit$format > 3){
        labels <- c('t','a0')
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
discordia_title <- function(fit,wetherill,sigdig=2,oerr=1,
                            y0option=1,dispunits=' Ma',...){
    if (is.null(fit$posterior) || 't'%ni%names(fit$posterior)){
        line1 <- maintit(x=fit$par[1],sx=fit$err[,1],n=fit$n,df=fit$df,
                         sigdig=sigdig,oerr=oerr,prefix='lower intercept =')
    } else {
        line1 <- bayestit(x=fit$par[1],XL=fit$posterior$t,n=fit$n,
                          sigdig=sigdig,oerr=oerr,prefix='lower intercept =')
    }
    if (wetherill){
        line2 <- maintit(x=fit$par[2],sx=fit$err[,2],ntit='',df=fit$df,
                         sigdig=sigdig,oerr=oerr,prefix='upper intercept =')
    } else if (fit$format %in% c(1,2,3)){
        if (is.null(fit$posterior)) pnames <- NULL
        else pnames <- names(fit$posterior)
        if (is.null(pnames)) ipar <- NULL
        else if (y0option==2 && 'U48i'%in%pnames) ipar <- 'U48i'
        else if (y0option==3 && 'ThUi'%in%pnames) ipar <-'ThUi'
        else ipar <- NULL
        fit <- getUPby0(fit,option=y0option)
        if (is.null(ipar)){
            line2 <- maintit(x=fit$par['a0'],sx=fit$err[,'a0'],ntit='',
                             sigdig=sigdig,oerr=oerr,units='',df=fit$df,
                             prefix=fit$y0label)
        } else {
            line2 <- bayestit(x=fit$par[ipar],XL=fit$posterior[[ipar]],ntit='',
                              sigdig=sigdig,oerr=oerr,units='',prefix=fit$y0label)
        }
    } else if (fit$format%in%c(4,5,6)){
        line2 <- maintit(x=fit$par['a0'],sx=fit$err[,'a0'],ntit='',
                         sigdig=sigdig,oerr=oerr,units='',df=fit$df,
                         prefix=quote('('^206*'Pb/'^204*'Pb)'[c]*'='))
        line3 <- maintit(x=fit$par['b0'],sx=fit$err[,'b0'],ntit='',
                         sigdig=sigdig,oerr=oerr,units='',df=fit$df,
                         prefix=quote('('^207*'Pb/'^204*'Pb)'[c]*'='))
    } else if (fit$format%in%c(7,8,85)){
        i86 <- 1/fit$par['a0']
        i87 <- 1/fit$par['b0']
        i86err <- i86*fit$err[,'a0']/fit$par['a0']
        i87err <- i87*fit$err[,'b0']/fit$par['b0']
        line2 <- maintit(x=i86,sx=i86err,ntit='',sigdig=sigdig,oerr=oerr,units='',
                         df=fit$df,prefix=quote('('^208*'Pb/'^206*'Pb)'[c]*'='))
        line3 <- maintit(x=i87,sx=i87err,ntit='',sigdig=sigdig,oerr=oerr,units='',
                         df=fit$df,prefix=quote('('^208*'Pb/'^207*'Pb)'[c]*'='))
    } else {
        stop('Invalid U-Pb data format.')
    }
    if (fit$model==1){
        line4 <- mswdtit(mswd=fit$mswd,p=fit$p.value,sigdig=sigdig)
    } else if (fit$model==3){
        line4 <- disptit(w=fit$disp['w'],sw=fit$disp['s[w]'],
                         units=dispunits,sigdig=sigdig,oerr=oerr)
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
