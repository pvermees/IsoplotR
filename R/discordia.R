# returns the lower and upper intercept age (for Wetherill concordia)
# or the lower intercept age and 207Pb/206Pb intercept (for Tera-Wasserburg)
concordia.intersection <- function(x,wetherill=TRUE,exterr=FALSE){
    if (x$format<4)
        out <- concordia.intersection.york(x,wetherill=wetherill,exterr=exterr)
    else
        out <- concordia.intersection.ludwig(x,wetherill=wetherill,exterr=exterr)
    out
}
concordia.intersection.york <- function(x,wetherill=TRUE,exterr=FALSE){
    d <- data2york(x,wetherill=wetherill)
    fit <- york(d)
    out <- list()
    if (wetherill){
        search.range <- c(0,10000)
        midpoint <- stats::optimize(intersection.misfit.york, search.range,
                                    a=fit$a[1], b=fit$b[1],
                                    wetherill=wetherill)$minimum
        range1 <- c(-1000,midpoint)
        range2 <- c(midpoint,10000)
        out$x <- search.range # tl, tu
        names(out$x) <- c('t[l]','t[u]')
        out$x['t[l]'] <- stats::uniroot(intersection.misfit.york, range1, 
                                        a=fit$a[1], b=fit$b[1], wetherill=wetherill)$root
        out$x['t[u]'] <- stats::uniroot(intersection.misfit.york, range2, 
                                        a=fit$a[1], b=fit$b[1], wetherill=wetherill)$root
    } else {
        search.range <- c(1/10000,10000)
        out$x <- c(1,fit$a[1]) # tl, 7/6 intercept
        names(out$x) <- c('t[l]','76')
        if (fit$b[1]<0) { # negative slope => two intersections with concordia line
            midpoint <- stats::optimize(intersection.misfit.york, search.range,
                                        a=fit$a[1], b=fit$b[1],
                                        wetherill=wetherill)$minimum
            search.range[2] <- midpoint
            out$x['t[l]'] <- stats::uniroot(intersection.misfit.york, search.range, 
                                            a=fit$a[1], b=fit$b[1],
                                            wetherill=wetherill)$root
        } else {
            out$x['t[l]'] <- stats::uniroot(intersection.misfit.york, search.range,
                                            a=fit$a[1], b=fit$b[1],
                                            wetherill=wetherill)$root
        }
    }
    hess <- stats::optimHess(out$x,LL.concordia.intersection.york,d=d,x=x,
                             wetherill=wetherill,exterr=exterr)
    out$cov <- solve(hess)
    out$mswd <- fit$mswd
    out$p.value <- fit$p.value
    out
}
# used by common Pb correction:
project.concordia <- function(m76,m86,i76){
    if (i76>m76)
        tend <- get.Pb206U238.age(1/m86)[1]
    else
        tend <- get.Pb207Pb206.age(m76)[1]
    search.range <- c(1/10000,tend)
    a <- i76
    b <- (m76-i76)/m86
    stats::uniroot(intersection.misfit.york, search.range, 
                   a=a,b=b,wetherill=FALSE)$root
}
concordia.intersection.ludwig <- function(x,wetherill=TRUE,exterr=FALSE){
    out <- list()
    fit <- ludwig(x,exterr=exterr)
    out$x <- c(0,0)
    J <- matrix(0,2,3)
    if (wetherill){
        t1 <- fit$par['t']
        a0 <- fit$par['64i']
        b0 <- fit$par['74i']
        names(out$x) <- c('t[l]','t[u]')
        buffer <- 1 # start searching 1Ma above or below first intercept age
        l5 <- lambda('U235')[1]
        l8 <- lambda('U238')[1]
        R <- iratio('U238U235')[1]
        disc.slope <- fit$par['64i']/(fit$par['74i']*R)
        conc.slope <- (l8*exp(l8*t1))/(l5*exp(l5*t1))
        if (conc.slope > disc.slope){
            search.range <- c(t1+buffer,10000)
            t1.name <- 't[l]'
            t2.name <- 't[u]'
        } else {
            search.range <- c(-1000,t1-buffer)
            t1.name <- 't[u]'
            t2.name <- 't[l]'
        }
        t2 <- stats::uniroot(intersection.misfit.ludwig,
                             interval=search.range,
                             t1=t1,a0=a0,b0=b0)$root
        out$x[t1.name] <- t1
        out$x[t2.name] <- t2
        J <- J.lud2york(t1,t2,a0,b0)
        out$cov <- J %*% fit$cov %*% t(J)
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
J.lud2york <- function(t1,t2,a0,b0){
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    R <- iratio('U238U235')[1]
    if (t2>t1){
        X <- exp(l5*t2) - exp(l5*t1)
        Y <- exp(l8*t2) - exp(l8*t1)
        dX.dt1 <- -l5*exp(l5*t1)
        dY.dt1 <- -l8*exp(l8*t1)
        dX.dt2 <- l5*exp(l8*t2)
        dY.dt2 <- l8*exp(l8*t2)
    } else {
        X <- exp(l5*t1) - exp(l5*t2)
        Y <- exp(l8*t1) - exp(l8*t2)
        dX.dt1 <- l5*exp(l5*t1)
        dY.dt1 <- l8*exp(l8*t1)
        dX.dt2 <- -l5*exp(l8*t2)
        dY.dt2 <- -l8*exp(l8*t2)
    }
    B <- a0/(b0*R)
    dB.da0 <- 1/(B*R)
    dB.db0 <- -B/b0
    dD.dt1 <- 2*(Y-B*X)*(dY.dt1-B*dX.dt1)
    dD.dt2 <- 2*(Y-B*X)*(dY.dt2-B*dX.dt2)
    dD.da0 <- 2*(B*X-Y)*dB.da0*X
    dD.db0 <- 2*(B*X-Y)*dB.db0*X
    J <- matrix(0,2,3)
    J[1,1] <- 1
    J[2,1] <- -dD.dt1/dD.dt2
    J[2,2] <- -dD.da0/dD.dt2
    J[2,3] <- -dD.db0/dD.dt2
    J
}

# itt = output of the york function
# d = output of UPb2york
# x = U-Pb data
LL.concordia.intersection.york <- function(itt,d,x,wetherill=TRUE,exterr=FALSE){
    LL <- 0
    if (wetherill){
        XYl <- age_to_wetherill_ratios(itt[1])
        XYu <- age_to_wetherill_ratios(itt[2])
        b <- (XYu$x[2]-XYl$x[2])/(XYu$x[1]-XYl$x[1])
        a <- XYu$x[2]-b*XYu$x[1]
    } else {
        XYl <- age_to_terawasserburg_ratios(itt[1])
        a <- itt[2]
        b <- (XYl$x[2]-a)/XYl$x[1]
    }
    xy <- get.york.xy(d,a,b)
    mu <- d[,c('X','Y')]-xy
    if (wetherill) disc <- (XYu$x[1]-xy[,1])/(XYu$x[1]-XYl$x[1])
    else disc <- (XYl$x[1]-xy[,1])/XYl$x[1]
    for (i in 1:length(x)){
        if (wetherill) samp <- wetherill(x,i)
        else samp <- tera.wasserburg(x,i)
        covmat <- samp$cov[1:2,1:2]
        if (exterr) {
            dcomp <- discordant.composition(disc[i],itt[1],itt[2],wetherill)
            covmat <- samp$cov[1:2,1:2] + dcomp$cov
        }
        LL <- LL + LL.norm(matrix(mu[i,],1,2),covmat)
    }
    LL
}
# returns upper intercept for for Wetherill concordia
LL.concordia.intersection.ludwig <- function(tu,tl,a0,b0){
    intersection.misfit.ludwig(tu,tl,a0,b0)^2
}

# returns misfit of a proposed age and the intersection between the
# discordia and concordia lines
intersection.misfit.york <- function(age,a,b,wetherill=TRUE){
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    R <- iratio('U238U235')[1]
    if (wetherill)
        out <- a-b+1 + b*exp(l5*age) - exp(l8*age)
    else
        out <- (exp(l5*age)-1)/(exp(l8*age)-1) - a*R - b*R/(exp(l8*age)-1)
    out
}
intersection.misfit.ludwig <- function(t2,t1,a0,b0){
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

# find the composition of a point that is 100d% discordant
# between upper intercept and lower intercept age
# (if wetherill=TRUE) or between the 207/206 intercept
# and the lower intercept age (if wetherill=FALSE)
discordant.composition <- function(d,tl,itu,wetherill=TRUE){
    out <- list()
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    R <- iratio('U238U235')[1]
    if (wetherill){
        X <- d*(exp(l5*tl)-1) + (1-d)*(exp(l5*itu)-1)
        Y <- d*(exp(l8*tl)-1) + (1-d)*(exp(l8*itu)-1)
        dXdl5 <- d*tl*exp(l5*tl) + (1-d)*itu*exp(l5*itu)
        dXdl8 <- 0
        dXdR <- 0
        dYdl5 <- 0
        dYdl8 <- d*tl*exp(l8*tl) + (1-d)*itu*exp(l8*itu)
        dYdR <- 0
    } else {
        X <- (1-d)/(exp(l8*tl)-1)
        Y <- d*itu + ((1-d)/R)*(exp(l5*tl)-1)/(exp(l8*tl)-1)
        dXdl5 <- 0
        dXdl8 <- -(1-d)*tl*exp(l8*tl)/(exp(l8*tl)-1)^2
        dXdR <- 0
        dYdl5 <- ((1-d)/R)*tl*exp(l5*tl)/(exp(l8*tl)-1)
        dYdl8 <- -((1-d)/R)*tl*exp(l8*tl)*(exp(l5*tl)-1)/(exp(l8*tl)-1)^2
        dYdR <- -((1-d)/R^2)*(exp(l5*tl)-1)/(exp(l8*tl)-1)
    }
    J <- matrix(0,2,3)
    J[1,1] <- dXdl5
    J[1,2] <- dXdl8
    J[1,3] <- dXdR
    J[2,1] <- dYdl5
    J[2,2] <- dYdl8
    J[2,3] <- dYdR
    E <- matrix(0,3,3)
    E[1,1] <- lambda('U235')[2]^2
    E[2,2] <- lambda('U238')[2]^2
    E[3,3] <- iratio('U238U235')[2]^2
    covmat <- J %*% E %*% t(J)
    out$x <- c(X,Y)
    out$cov <- covmat
    out
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
