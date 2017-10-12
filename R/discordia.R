# returns the lower and upper intercept age (for Wetherill concordia)
# or the lower intercept age and 207Pb/206Pb intercept (for Tera-Wasserburg)
concordia.intersection.ludwig <- function(x,wetherill=TRUE,
                                          exterr=FALSE,alpha=0.05,
                                          model=1){
    fit <- ludwig(x,exterr=exterr,model=model)
    out <- list()
    out$model <- model
    out$mswd <- fit$mswd
    out$w <- fit$w
    out$p.value <- fit$p.value
    out$df <- fit$df
    tfact <- stats::qt(1-alpha/2,fit$df)
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
        out$err['disp',] <- tfact*sqrt(fit$mswd)*sqrt(diag(out$cov))
    } else {
        out$err <- matrix(NA,2,2)
        rownames(out$err) <- c('s','ci')
    }
    out$err['s',] <- sqrt(diag(out$cov))
    out$err['ci',] <- tfact*out$err['s',]
    colnames(out$err) <- labels
    out
}
concordia.intersection.york <- function(x,exterr=FALSE){
    d <- data2york(x,wetherill=FALSE)
    fit <- york(d)
    concordia.intersection.ab(fit$a[1],fit$b[1],exterr=exterr)
}
concordia.intersection.ab <- function(a,b,exterr=FALSE,wetherill=FALSE){
    out <- list()
    m <- 1/10000
    M <- 10000
    out$x <- c(1,a) # tl, 7/6 intercept
    if (wetherill) names(out$x) <- c('t[l]','t[u]')
    else names(out$x) <- c('t[l]','76i')
    if (b<0) { # negative slope => two intersections with concordia line
        search.range <- c(m,M)
        midpoint <- stats::optimize(intersection.misfit.york,
                                    search.range,a=a,b=b)$minimum
        search.range <- c(m,midpoint)
        out$x['t[l]'] <- stats::uniroot(intersection.misfit.york,
                                        search.range,a=a,b=b)$root
        if (wetherill){
            search.range <- c(midpoint,M)
            out$x['t[u]'] <- stats::uniroot(intersection.misfit.york,
                                            search.range,a=a,b=b)$root
        }   
    } else {
        search.range <- c(m,M)
        out$x['t[l]'] <- stats::uniroot(intersection.misfit.york,
                                        search.range,a=a,b=b)$root
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
    R <- iratio('U238U235')[1]
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
    disc.slope <- a0/(b0*R)
    conc.slope <- (l8*exp(l8*tt))/(l5*exp(l5*tt))
    if (conc.slope > disc.slope){
        search.range <- c(tt+buffer,10000)
        tl <- tt
        tu <- stats::uniroot(intersection.misfit.ludwig,
                             interval=search.range,
                             t2=tt,a0=a0,b0=b0)$root
    } else {
        search.range <- c(-1000,tt-buffer)
        tl <- stats::uniroot(intersection.misfit.ludwig,
                             interval=search.range,
                             t2=tt,a0=a0,b0=b0)$root
        tu <- tt
    }
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    R <- iratio('U238U235')[1]
    XX <- exp(l5*tu) - exp(l5*tl)
    YY <- exp(l8*tu) - exp(l8*tl)
    BB <- a0/(b0*R)
    D <- (YY-BB*XX)^2
    dXX.dtl <- -l5*exp(l5*tl)
    dXX.dtu <-  l5*exp(l5*tu)
    dYY.dtl <- -l8*exp(l8*tl)
    dYY.dtu <-  l8*exp(l8*tu)
    dBB.da0 <-  1/(b0*R)
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

# used by common Pb correction:
project.concordia <- function(m76,m86,i76){
    t68 <- get.Pb206U238.age(1/m86)[1]
    t76 <- get.Pb207Pb206.age(m76)[1]
    a <- i76
    b <- (m76-i76)/m86
    neg <- (i76>m76) # negative slope?
    pos <- !neg
    above <- (t76>t68) # above concordia?
    below <- !above
    tend <- t68
    go.ahead <- FALSE
    if (pos & above){
        tend <- t76
        go.ahead <- TRUE
    } else if (pos & below){
        go.ahead <- TRUE
    } else if (neg & above){
        go.ahead <- TRUE
    } else if (neg & below){  # it is not clear what to do with samples
        for (tt in seq(from=10,to=5000,by=10)){
            misfit <- intersection.misfit.york(tt,a=a,b=b)
            if (misfit<0){    # that plot in the 'forbidden zone' above
                tend <- tt    # Wetherill concordia or below T-W concordia
                go.ahead <- TRUE # IsoplotR will still project them on
                break         # the concordia line.
            }
        }
    }
    if (go.ahead){
        search.range <- c(1/10000,tend)
        out <- stats::uniroot(intersection.misfit.york,search.range,
                              a=a,b=b)$root
    } else {
        out <- m76
    }
    out
}

# returns misfit of a proposed age and the intersection
# between the discordia and concordia lines
intersection.misfit.ludwig <- function(t1,t2,a0,b0){
    tl <- min(t1,t2)
    tu <- max(t1,t2)
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    R <- iratio('U238U235')[1]
    XX <- exp(l5*tu) - exp(l5*tl)
    YY <- exp(l8*tu) - exp(l8*tl)
    BB <- a0/(b0*R)
    YY - BB*XX
}
intersection.misfit.york <- function(age,a,b){
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    R <- iratio('U238U235')[1]
    (exp(l5*age)-1)/(exp(l8*age)-1) - a*R - b*R/(exp(l8*age)-1)
}

discordia.line <- function(fit,wetherill){
    X <- c(0,0)
    Y <- c(0,0)
    if (wetherill){
        X <- age_to_Pb207U235_ratio(fit$x)[,'75']
        Y <- age_to_Pb206U238_ratio(fit$x)[,'68']
    } else {
        X[1] <- age_to_U238Pb206_ratio(fit$x['t[l]'])[,'86']
        Y[1] <- age_to_Pb207Pb206_ratio(fit$x['t[l]'])[,'76']
        Y[2] <- fit$x['76']
    }
    graphics::lines(X,Y)
}

# this would be much easier in unicode but that doesn't render in PDF:
discordia.title <- function(fit,wetherill,sigdig=2){
    lower.age <- roundit(fit$x[1],fit$err[,1],sigdig=sigdig)
    list1 <- list(a=lower.age[1],b=lower.age[2],c=lower.age[3])
    if (fit$model!=2 && fit$mswd>1) args <- quote(a%+-%b~'|'~c~'|'~d)
    else args <- quote(a%+-%b~'|'~c)
    if (wetherill){
        upper.age <- roundit(fit$x[2],fit$err[,2],sigdig=sigdig)
        expr1 <- quote('lower intercept =')
        expr2 <- quote('upper intercept =')
        list2 <- list(a=upper.age[1],b=upper.age[2],c=upper.age[3])
        if (fit$model!=2 && fit$mswd>1){
            list1$d <- lower.age[4]
            list2$d <- upper.age[4]
        }
    } else {
        intercept <- roundit(fit$x[2],fit$err[,2],sigdig=sigdig)
        expr1 <- quote('age =')
        expr2 <- quote('('^207*'Pb/'^206*'Pb)'[o]~'=')
        list2 <- list(a=intercept[1],b=intercept[2],c=intercept[3])
        if (fit$model!=2 && fit$mswd>1){
            list1$d <- lower.age[4]
            list2$d <- intercept[4]
        }
    }
    call1 <- substitute(e~a,list(e=expr1,a=args))
    call2 <- substitute(e~a,list(e=expr2,a=args))
    line1 <- do.call('substitute',list((call1),list1))
    line2 <- do.call('substitute',list((call2),list2))
    if (fit$model==1){
        line3 <- substitute('MSWD ='~a~', p('~chi^2*')='~b,
                            list(a=signif(fit$mswd,sigdig),
                                 b=signif(fit$p.value,sigdig)))
        graphics::mtext(line1,line=2)
        graphics::mtext(line2,line=1)
        graphics::mtext(line3,line=0)
    } else if (fit$model==2){
        graphics::mtext(line1,line=1)
        graphics::mtext(line2,line=0)        
    } else {
        line3 <- substitute('overdispersion ='~a~'|'~b~
                            '% of the initial Pb',
                            list(a=signif(100*fit$w['s'],sigdig),
                                 b=signif(100*fit$w['ci'],sigdig)))
        graphics::mtext(line1,line=2)
        graphics::mtext(line2,line=1)
        graphics::mtext(line3,line=0)
    }
}
