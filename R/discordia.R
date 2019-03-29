# returns the lower and upper intercept age (for Wetherill concordia)
# or the lower intercept age and 207Pb/206Pb intercept (for Tera-Wasserburg)
concordia.intersection.ludwig <- function(x,wetherill=TRUE,exterr=FALSE,alpha=0.05,
                                          model=1,anchor=list(FALSE,NA)){
    fit <- ludwig(x,exterr=exterr,model=model,anchor=anchor)
    out <- list()
    out$model <- model
    out$mswd <- fit$mswd
    out$w <- fit$w
    out$p.value <- fit$p.value
    out$fact <- tfact(alpha,fit$df)
    if (wetherill){
        labels <- c('t[l]','t[u]')
        out <- c(out,twfit2wfit(fit,x))
    } else if (x$format<4){
        labels <- c('t[l]','76')
        out$x <- fit$par
        out$cov <- fit$cov
    } else {
        labels <- c('t[l]','76')
        out$x <- c(fit$par['t'],
                   fit$par['74i']/fit$par['64i'])
        J <- matrix(0,2,3)
        J[1,1] <- 1
        J[2,2] <- -fit$par['74i']/fit$par['64i']^2
        J[2,3] <- 1/fit$par['64i']
        out$cov <- J %*% fit$cov %*% t(J)
    }
    names(out$x) <- labels
    if (model==1 && fit$mswd>1){
        out$err <- matrix(NA,3,2)
        rownames(out$err) <- c('s','ci','disp')
        out$err['disp',] <-
            out$fact*sqrt(fit$mswd)*sqrt(diag(out$cov))
    } else {
        out$err <- matrix(NA,2,2)
        rownames(out$err) <- c('s','ci')
    }
    if (model==3) out$fact <- nfact(alpha)
    out$err['s',] <- sqrt(diag(out$cov))
    out$err['ci',] <- out$fact*out$err['s',]
    colnames(out$err) <- labels
    out
}
# extracts concordia intersection parameters from an ordinary York fit
concordia.intersection.ab <- function(a,b,exterr=FALSE,wetherill=FALSE,d=diseq()){
    l8 <- lambda('U238')[1]
    ta <- get.Pb207Pb206.age(a,d=d)[1]
    out <- c(1,a) # tl, 7/6 intercept
    if (wetherill) names(out) <- c('t[l]','t[u]')
    else names(out) <- c('t[l]','76i')
    if (b<0) { # negative slope => two intersections with concordia line
        tb <- get.Pb206U238.age(-b/a,d=d)[1]
        search.range <- c(tb,ta)
        midpoint <- stats::optimize(intersection.misfit.york,
                                    search.range,a=a,b=b,d=d)$minimum
        search.range[2] <- midpoint
        out['t[l]'] <- stats::uniroot(intersection.misfit.york,
                                      interval=search.range,a=a,b=b,d=d)$root
        if (wetherill){
            search.range <- c(midpoint,ta)
            out['t[u]'] <- stats::uniroot(intersection.misfit.york,
                                          interval=search.range,a=a,b=b,d=d)$root
        }   
    } else {
        search.range <- c(ta,2/l8)
        out['t[l]'] <- stats::uniroot(intersection.misfit.york,
                                      interval=search.range,a=a,b=b,d=d)$root
    }
    out
}

# extract the lower and upper discordia intercept from the parameters
# of a Ludwig fit (initial Pb ratio and lower intercept age)
twfit2wfit <- function(fit,x){
    tt <- fit$par[1]
    buffer <- 1 # start searching 1Ma above or below first intercept age
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    U <- iratio('U238U235')[1]
    E <- matrix(0,3,3)
    J <- matrix(0,2,3)
    if (x$format<4){
        a0 <- 1
        b0 <- fit$par['76i']
        E[c(1,3),c(1,3)] <- fit$cov
    } else {
        a0 <- fit$par['64i']
        b0 <- fit$par['74i']
        E <- fit$cov
    }
    disc.slope <- a0/(b0*U)
    conc.slope <- (l8*exp(l8*tt))/(l5*exp(l5*tt))
    if (conc.slope > disc.slope){
        search.range <- c(tt,get.Pb207Pb206.age(b0/a0,d=x$d)[1])+buffer
        tl <- tt
        tu <- stats::uniroot(intersection.misfit.ludwig,interval=search.range,
                             t2=tt,a0=a0,b0=b0,d=x$d)$root
    } else {
        search.range <- c(-1000,tt-buffer)
        if (x$d$corr) search.range[1] <- -0.001/settings('lambda','U234')[1]
        tl <- stats::uniroot(intersection.misfit.ludwig,interval=search.range,
                             t2=tt,a0=a0,b0=b0,d=x$d)$root
        tu <- tt
    }
    du <- wendt(tt=tu,d=x$d)
    dl <- wendt(tt=tl,d=x$d)
    XX <- exp(l5*tu) + du$d1 - exp(l5*tl) - dl$d1
    YY <- exp(l8*tu) + du$d2 - exp(l8*tl) - dl$d2
    BB <- a0/(b0*U)
    D <- (YY-BB*XX)^2 # misfit
    dXX.dtl <- -l5*exp(l5*tl) - dl$dd1dt
    dXX.dtu <-  l5*exp(l5*tu) + du$dd1dt
    dYY.dtl <- -l8*exp(l8*tl) - dl$dd2dt
    dYY.dtu <-  l8*exp(l8*tu) + du$dd2dt
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
    out$x <- c(tl,tu)
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
    du <- wendt(tt=tu,d=d)
    dl <- wendt(tt=tl,d=d)
    XX <- exp(l5*tu) + du$d1 - exp(l5*tl) - dl$d1
    YY <- exp(l8*tu) + du$d2 - exp(l8*tl) - dl$d2
    BB <- a0/(b0*U)
    # misfit is based on difference in slope in Wetherill space
    YY - BB*XX
}
# a = intercept, b = slope on TW concordia
intersection.misfit.york <- function(tt,a,b,d=diseq()){
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    U <- iratio('U238U235')[1]
    D <- wendt(tt=tt,d=d)
    # misfit is based on difference in slope in TW space
    (exp(l5*tt)-1+D$d1)/U - a*(exp(l8*tt)-1+D$d2) - b
}

discordia.line <- function(fit,wetherill,d=diseq()){
    X <- c(0,0)
    Y <- c(0,0)
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    J <- matrix(0,1,2)
    usr <- graphics::par('usr')
    if (wetherill){
        tl <- fit$x[1]
        tu <- fit$x[2]
        X <- age_to_Pb207U235_ratio(fit$x,d=d)[,'75']
        Y <- age_to_Pb206U238_ratio(fit$x,d=d)[,'68']
        x <- seq(from=max(0,usr[1],X[1]),to=min(usr[2],X[2]),length.out=50)
        du <- wendt(tt=tu,d=d)
        dl <- wendt(tt=tl,d=d)
        aa <- exp(l8*tu) + du$d2 - exp(l8*tl) - dl$d2
        bb <- x - exp(l5*tl) + 1 - dl$d1
        cc <- exp(l5*tu) + du$d1 - exp(l5*tl) - dl$d1
        dd <- exp(l8*tl) - 1 + dl$d2
        y <- aa*bb/cc + dd
        dadtl <- - l8*exp(l8*tl) - dl$dd2dt
        dbdtl <- - l5*exp(l5*tl) - dl$dd1dt
        dcdtl <- - l5*exp(l5*tl) - dl$dd1dt
        dddtl <- l8*exp(l8*tl) + dl$dd2dt
        dadtu <- l8*exp(l8*tu) + du$dd2dt
        dbdtu <- 0
        dcdtu <- l5*exp(l5*tu) + du$dd1dt
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
    } else {
        X[1] <- age_to_U238Pb206_ratio(fit$x['t[l]'],d=d)[,'86']
        Y[1] <- age_to_Pb207Pb206_ratio(fit$x['t[l]'],d=d)[,'76']
        r75 <- age_to_Pb207U235_ratio(fit$x['t[l]'],d=d)[,'75']
        r68 <- 1/X[1]
        Y[2] <- fit$x['76']
        xl <- X[1]
        yl <- Y[1]
        y0 <- Y[2]
        tl <- check.zero.UPb(fit$x['t[l]'])
        U <- settings('iratio','U238U235')[1]
        nsteps <- 100
        x <- seq(from=max(.Machine$double.xmin,usr[1]),to=usr[2],length.out=nsteps)
        y <- yl + (y0-yl)*(1-x*r68) # = y0 + yl*x*r68 - y0*x*r68
        D <- wendt(tt=tl,d=d)
        d75dtl <- l5*exp(l5*tl) + D$dd1dt
        d68dtl <- l8*exp(l8*tl) + D$dd2dt
        dyldtl <- (d75dtl*r68 - r75*d68dtl)/(U*r68^2)
        J1 <- dyldtl*x*r68 + yl*x*d68dtl - y0*x*d68dtl # dy/dtl
        J2 <- 1 - x*r68                                # dy/dy0
        sy <- errorprop1x2(J1,J2,fit$cov[1,1],fit$cov[2,2],fit$cov[1,2])
        ul <- y + fit$fact*sy
        ll <- y - fit$fact*sy
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

# this would be much easier in unicode but that doesn't render in PDF:
discordia.title <- function(fit,wetherill,sigdig=2,...){
    lower.age <- roundit(fit$x[1],fit$err[,1],sigdig=sigdig)
    if (fit$model==1 && fit$mswd>1){
        args1 <- quote(a%+-%b~'|'~c~'|'~d~u~'(n='*n*')')
        args2 <- quote(a%+-%b~'|'~c~'|'~d~u)
    } else {
        args1 <- quote(a%+-%b~'|'~c~u~'(n='*n*')')
        args2 <- quote(a%+-%b~'|'~c~u)
    }
    list1 <- list(a=lower.age[1],b=lower.age[2],
                  c=lower.age[3],u='Ma',n=fit$n)
    if (wetherill){
        upper.age <- roundit(fit$x[2],fit$err[,2],sigdig=sigdig)
        expr1 <- quote('lower intercept =')
        expr2 <- quote('upper intercept =')
        list2 <- list(a=upper.age[1],b=upper.age[2],c=upper.age[3],u='Ma')
        if (fit$model==1 && fit$mswd>1){
            list1$d <- lower.age[4]
            list2$d <- upper.age[4]
        }
    } else {
        intercept <- roundit(fit$x[2],fit$err[,2],sigdig=sigdig)
        expr1 <- quote('age =')
        expr2 <- quote('('^207*'Pb/'^206*'Pb)'[o]~'=')
        list2 <- list(a=intercept[1],b=intercept[2],c=intercept[3],u='')
        if (fit$model==1 && fit$mswd>1){
            list1$d <- lower.age[4]
            list2$d <- intercept[4]
        }
    }
    call1 <- substitute(e~a,list(e=expr1,a=args1))
    call2 <- substitute(e~a,list(e=expr2,a=args2))
    line1 <- do.call('substitute',list((call1),list1))
    line2 <- do.call('substitute',list((call2),list2))
    if (fit$model==1){
        line3 <- substitute('MSWD ='~a~', p('*chi^2*') ='~b,
                            list(a=signif(fit$mswd,sigdig),
                                 b=signif(fit$p.value,sigdig)))
        mymtext(line1,line=2,...)
        mymtext(line2,line=1,...)
        mymtext(line3,line=0,...)
    } else if (fit$model==2){
        mymtext(line1,line=1,...)
        mymtext(line2,line=0,...)
    } else {
        rounded.disp <- roundit(100*fit$w[1],100*fit$w[2:3],sigdig=sigdig)
        line3 <- substitute('overdispersion ='~a+b/-c~
                            '% of Pb'[o],
                            list(a=rounded.disp[1],
                                 b=rounded.disp[3],
                                 c=rounded.disp[2]))
        mymtext(line1,line=2,...)
        mymtext(line2,line=1,...)
        mymtext(line3,line=0,...)
    }
}
