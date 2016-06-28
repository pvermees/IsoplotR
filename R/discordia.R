discordia.age <- function(x,wetherill=TRUE,dcu=TRUE){
    out <- list()
    selection <- get.selection(x,wetherill)
    d <- data2york(x,selection)
    fit <- yorkfit(d$X,d$Y,d$sX,d$sY,d$rXY)
    itt <- concordia.intersection(fit,wetherill)
    X <- UPb.preprocess(x,wetherill)
    hess <- stats::optimHess(itt,LL.concordia.intersection, d=d,X=X,
                             wetherill=wetherill,dcu=dcu)
    out$x <- itt
    out$cov <- solve(hess)
    out
}

# itt = output of the yorkfit function
# d = output of UPb2york
# X = output of UPb.preprocess
LL.concordia.intersection <- function(itt,d,X,wetherill,dcu){
    LL <- 0
    selection <- get.UPb.labels(wetherill)
    XYl <- get.ratios.UPb(itt[1])$x[selection]
    if (wetherill){
        XYu <- get.ratios.UPb(itt[2])$x[selection]
        b <- (XYu[2]-XYl[2])/(XYu[1]-XYl[1])
        a <- XYu[2]-b*XYu[1]
    } else {
        a <- itt[2]
        b <- (XYl[2]-a)/XYl[1]
    }
    xy <- get.york.xy(d$X,d$Y,d$sX,d$sY,d$rXY,a,b)
    mu <- cbind(d$X,d$Y)-xy
    if (wetherill) disc <- (XYu[1]-xy[,1])/(XYu[1]-XYl[1])
    else disc <- (XYl[1]-xy[,1])/XYl[1]
    for (i in 1:length(X)){
        covmat <- X[[i]]$cov
        if (dcu) {
            dcomp <- discordant.composition(disc[i],itt[1],itt[2],wetherill)
            covmat <- covmat + dcomp$cov
        }
        LL <- LL + LL.norm(matrix(mu[i,],1,2),covmat)
    }
    LL
}

# returns the lower and upper intercept age (for Wetherill concordia)
# or the lower intercept age and 207Pb/206Pb intercept (for Tera-Wasserburg)
concordia.intersection <- function(fit,wetherill){
    if (wetherill){
        search.range <- c(0,10000)
        midpoint <- stats::optimize(intersection.misfit, search.range,
                                    a=fit$a[1], b=fit$b[1], wetherill=wetherill)$minimum
        range1 <- c(-1000,midpoint)
        range2 <- c(midpoint,10000)
        tt <- search.range # tl, tu
        tt[1] <- stats::uniroot(intersection.misfit, range1, 
                                a=fit$a[1], b=fit$b[1], wetherill=wetherill)$root
        tt[2] <- stats::uniroot(intersection.misfit, range2, 
                                a=fit$a[1], b=fit$b[1], wetherill=wetherill)$root
        out <- tt
    } else {
        search.range <- c(1/10000,10000)
        it <- c(1,fit$a[1]) # tl, 7/6 intercept
        if (fit$b[1]<0) { # negative slope => two intersections with concordia line
            midpoint <- stats::optimize(intersection.misfit, search.range,
                                        a=fit$a[1], b=fit$b[1], wetherill=wetherill)$minimum
            search.range[2] <- midpoint
            it[1] <- stats::uniroot(intersection.misfit, search.range, 
                                    a=fit$a[1], b=fit$b[1], wetherill=wetherill)$root
            out <- it
        } else {
            it[1] <- stats::uniroot(intersection.misfit, search.range,
                                    a=fit$a[1], b=fit$b[1], wetherill=wetherill)$root
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
    if (wetherill){
        out <- a-b+1 + b*exp(l5*age) - exp(l8*age)
    } else {
        out <- (exp(l5*age)-1)/(exp(l8*age)-1) - a*R - b*R/(exp(l8*age)-1)
    }
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
    labels <- get.UPb.labels(wetherill)
    names(out$x) <- labels
    rownames(out$cov) <- labels
    colnames(out$cov) <- labels
    out
}

# negative multivariate log likelihood to be fed into R's optim function
LL.norm <- function(x,covmat){
    log(2*pi) + 0.5*determinant(covmat,logarithmic=TRUE)$modulus + 0.5*get.SS(x,covmat)
}

get.SS <- function(x,covmat){
    x %*% solve(covmat) %*% t(x)
}

discordia.plot <- function(fit,wetherill=TRUE){
    selection <- get.UPb.labels(wetherill)
    comp.l <- get.ratios.UPb(fit$x[1])$x[selection]
    if (wetherill) {
        comp.u <- get.ratios.UPb(fit$x[2])$x[selection]
    } else {
        comp.u <- c(0,fit$x[2])
    }
    X <- c(comp.l[1],comp.u[1])
    Y <- c(comp.l[2],comp.u[2])
    graphics::lines(X,Y)
}

discordia.title <- function(fit,wetherill){
    if (wetherill){
        lower.age <- roundit(fit$x[1],sqrt(fit$cov[1,1]))
        upper.age <- roundit(fit$x[2],sqrt(fit$cov[2,2]))
        line1 <- substitute('lower intercept ='~a%+-%b~'[Ma]',
                            list(a=lower.age$x, b=lower.age$err))
        line2 <- substitute('upper intercept ='~a%+-%b~'[Ma]',
                            list(a=upper.age$x, b=upper.age$err))
    } else {
        lower.age <- roundit(fit$x[1],sqrt(fit$cov[1,1]))
        intercept <- roundit(fit$x[2],sqrt(fit$cov[2,2]))
        line1 <- substitute('age ='~a%+-%b~'[Ma]',
                            list(a=lower.age$x, b=lower.age$err))
        line2 <- substitute('('^207*'Pb/'^206*'Pb)'[0]~'='~a%+-%b,
                              list(a=intercept$x, b=intercept$err))
    }
    graphics::mtext(line1,line=1)
    graphics::mtext(line2,line=0)
}
