# returns the lower and upper intercept age (for Wetherill concordia)
# or the lower intercept age and 207Pb/206Pb intercept (for Tera-Wasserburg)
concordia.intersection.ludwig <- function(x,wetherill=TRUE,exterr=FALSE){
    out <- list()
    fit <- ludwig(x,exterr=exterr)
    out$x <- c(0,0)
    J <- matrix(0,2,3)
    E <- matrix(0,3,3)
    if (wetherill){
        tt <- fit$par[1]
        names(out$x) <- c('t[l]','t[u]')
        buffer <- 1 # start searching 1Ma above or below first intercept age
        l5 <- lambda('U235')[1]
        l8 <- lambda('U238')[1]
        R <- iratio('U238U235')[1]
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
        out$x['t[l]'] <- tl
        out$x['t[u]'] <- tu
        out$cov <- J %*% E %*% t(J)
    } else if (x$format<4){
        out$x <- fit$par
        names(out$x) <- c('t[l]','76')
        out$cov <- fit$cov
    } else {
        names(out$x) <- c('t[l]','76')
        out$x['t[l]'] <- fit$par['t']
        out$x['76'] <- fit$par['74i']/fit$par['64i']
        J[1,1] <- 1
        J[2,2] <- -out$x['76']/fit$par['64i']
        J[2,3] <- 1/fit$par['64i']
        out$cov <- J %*% fit$cov %*% t(J)
    }
    out$mswd <- fit$mswd
    out$p.value <- fit$p.value
    out
}
concordia.intersection.york <- function(x,exterr=FALSE){
    d <- data2york(x,wetherill=FALSE)
    fit <- york(d)
    concordia.intersection.york.ab(fit$a[1],fit$b[1],exterr=exterr)
}
concordia.intersection.york.ab <- function(a,b,exterr=FALSE){
    out <- list()
    search.range <- c(1/10000,10000)
    out$x <- c(1,a) # tl, 7/6 intercept
    names(out$x) <- c('t[l]','76i')
    if (b<0) { # negative slope => two intersections with concordia line
        midpoint <- stats::optimize(intersection.misfit.york,
                                    search.range,a=a,b=b)$minimum
        search.range[2] <- midpoint
        out$x['t[l]'] <- stats::uniroot(intersection.misfit.york,
                                        search.range,a=a,b=b)$root
    } else {
        out$x['t[l]'] <- stats::uniroot(intersection.misfit.york,
                                        search.range,a=a,b=b)$root
    }
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
            misfit <- intersection.misfit.york(tend,a=a,b=b)
            if (misfit<0){    # that plot in the 'forbidden zone' above
                tend <- tt    # Wetherill concordia or below T-W concordia
                found <- TRUE # IsoplotR will still project them on
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

discordia.plot <- function(fit,wetherill){
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

discordia.title <- function(fit,wetherill,sigdig=2){
    if (wetherill){
        lower.age <- roundit(fit$x[1],sqrt(fit$cov[1,1]),sigdig=sigdig)
        upper.age <- roundit(fit$x[2],sqrt(fit$cov[2,2]),sigdig=sigdig)
        line1 <- substitute('lower intercept ='~a%+-%b~'[Ma]',
                            list(a=lower.age[1], b=lower.age[2]))
        line2 <- substitute('upper intercept ='~a%+-%b~'[Ma]',
                            list(a=upper.age[1], b=upper.age[2]))
    } else {
        lower.age <- roundit(fit$x[1],sqrt(fit$cov[1,1]),sigdig=sigdig)
        intercept <- roundit(fit$x[2],sqrt(fit$cov[2,2]),sigdig=sigdig)
        line1 <- substitute('age ='~a%+-%b~'[Ma]',
                            list(a=lower.age[1], b=lower.age[2]))
        line2 <- substitute('('^207*'Pb/'^206*'Pb)'[o]~'='~a%+-%b,
                            list(a=intercept[1], b=intercept[2]))
    }
    line3 <- substitute('MSWD ='~a~', p('~chi^2*')='~b,
                        list(a=signif(fit$mswd,sigdig),
                             b=signif(fit$p.value,sigdig)))
    graphics::mtext(line1,line=2)
    graphics::mtext(line2,line=1)
    graphics::mtext(line3,line=0)
}
