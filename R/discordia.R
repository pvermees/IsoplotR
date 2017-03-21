discordia.age <- function(x,wetherill=TRUE,exterr=TRUE){
    d <- data2york(x,wetherill=wetherill)
    fit <- yorkfit(d)
    itt <- concordia.intersection(fit,wetherill)
    hess <- stats::optimHess(itt,LL.concordia.intersection, d=d,x=x,
                             wetherill=wetherill,exterr=exterr)
    out <- list()
    out$x <- itt
    out$cov <- solve(hess)
    out
}

# itt = output of the yorkfit function
# d = output of UPb2york
# X = output of UPb.preprocess
LL.concordia.intersection <- function(itt,d,x,wetherill,exterr){
    LL <- 0
    if (wetherill){
        XYl <- age.to.wetherill.ratios(itt[1])
        XYu <- age.to.wetherill.ratios(itt[2])
        b <- (XYu$x[2]-XYl$x[2])/(XYu$x[1]-XYl$x[1])
        a <- XYu$x[2]-b*XYu$x[1]
    } else {
        XYl <- age.to.terawasserburg.ratios(itt[1])
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
        if (exterr) {
            dcomp <- discordant.composition(disc[i],itt[1],itt[2],wetherill)
            covmat <- samp$cov + dcomp$cov
        }
        LL <- LL + LL.norm(matrix(mu[i,],1,2),covmat)
    }
    LL
}

# returns the lower and upper intercept age (for Wetherill concordia)
# or the lower intercept age and 207Pb/206Pb intercept (for Tera-Wasserburg)
concordia.intersection <- function(fit,wetherill=TRUE){
    if (wetherill){
        search.range <- c(0,10000)
        midpoint <- stats::optimize(intersection.misfit, search.range,
                                    a=fit$a[1], b=fit$b[1],
                                    wetherill=wetherill)$minimum
        range1 <- c(-1000,midpoint)
        range2 <- c(midpoint,10000)
        tt <- search.range # tl, tu
        tt[1] <- stats::uniroot(intersection.misfit, range1, 
                                a=fit$a[1], b=fit$b[1], wetherill=wetherill)$root
        tt[2] <- stats::uniroot(intersection.misfit, range2, 
                                a=fit$a[1], b=fit$b[1], wetherill=wetherill)$root
        out <- tt
        names(out) <- c('t[l]','t[u]')
    } else {
        search.range <- c(1/10000,10000)
        it <- c(1,fit$a[1]) # tl, 7/6 intercept
        names(it) <- c('t[l]','Pb207Pb206')
        if (fit$b[1]<0) { # negative slope => two intersections with concordia line
            midpoint <- stats::optimize(intersection.misfit, search.range,
                                        a=fit$a[1], b=fit$b[1],
                                        wetherill=wetherill)$minimum
            search.range[2] <- midpoint
            it[1] <- stats::uniroot(intersection.misfit, search.range, 
                                    a=fit$a[1], b=fit$b[1],
                                    wetherill=wetherill)$root
        } else {
            it[1] <- stats::uniroot(intersection.misfit, search.range,
                                    a=fit$a[1], b=fit$b[1],
                                    wetherill=wetherill)$root
        }
        out <- it
    }
    out
}

# returns misfit of a proposed age and the intersection between the
# discordia and concordia lines
intersection.misfit <- function(age,a,b,wetherill){
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    R <- iratio('U238U235')[1]
    if (wetherill)
        out <- a-b+1 + b*exp(l5*age) - exp(l8*age)
    else
        out <- (exp(l5*age)-1)/(exp(l8*age)-1) - a*R - b*R/(exp(l8*age)-1)
    out
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

discordia.plot <- function(fit,wetherill=TRUE){
    X <- c(0,0)
    Y <- c(0,0)
    if (wetherill) {
        X <- age.to.Pb207U235.ratio(fit$x)[,'75']
        Y <- age.to.Pb206U238.ratio(fit$x)[,'68']
    } else {
        X[1] <- age.to.U238Pb206.ratio(fit$x['t[l]'])[,'86']
        Y[1] <- age.to.Pb207Pb206.ratio(fit$x['t[l]'])[,'76']
        Y[2] <- fit$x['Pb207Pb206']
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
        line2 <- substitute('('^207*'Pb/'^206*'Pb)'[0]~'='~a%+-%b,
                              list(a=intercept[1], b=intercept[2]))
    }
    graphics::mtext(line1,line=1)
    graphics::mtext(line2,line=0)
}
