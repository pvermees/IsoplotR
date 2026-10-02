#' Leftmost and rightmost isochrons
#'
#' Modifies the Galbraith and Laslett's minimum age module to fit
#' overdispersed Pb-Pb and Th-U isochron data.
#' @param x an IsoplotR data object
#' @param left logical, switches between leftmost and rightmost
#'     isochrons
#' @param hide vector with indices of aliquots that should be removed
#'     from the plot.
#' @param omit vector with indices of aliquots that should be plotted
#'     but omitted from the isochron age calculation.
#' @param inverse toggles between normal and inverse isochrons. See
#'     \code{\link{isochron}}.
#' @param ... unused optional arguments
#' @return a list with the intercept (\code{a}), slope (\code{b}) and
#'     other data that can be passed on to \code{scatterplot}.
#' @noRd
LRisochron <- function(x,...){ UseMethod("LRisochron",x) }
#' @noRd
LRisochron.default <- function(x,left=TRUE,
                               hide=NULL,omit=NULL,
                               a=NULL,b=NULL,...){
    x2calc <- clear(x,hide,omit)
    sgn <- ifelse(left,-1,1)
    x2calc[,'Y'] <- sgn*x2calc[,'Y']
    X <- x2calc[,'X']
    Y <- x2calc[,'Y']
    h <- chull(X,Y)
    num_vertices <- length(h)
    hull_indices <- c(h,h[1])
    anchor <- 0
    if (is.null(b)){
        b <- diff(Y[hull_indices])/diff(X[hull_indices])
    } else {
        anchor <- 'b'
        b <- rep(b,num_vertices)
    }
    if (is.null(a)){
        a <- Y[h] - b*X[h]
    } else {
        anchor <- 'a'
        a <- rep(a,num_vertices)
    }
    lpi <- 0
    misfit <- Inf
    bestfit <- NULL
    for (i in 1:num_vertices){
        lsi <- log(sd((Y-a[i])/X))
        init <- c(lpi,lsi,a[i],b[i])
        testfit <- stats::optim(par=init,fn=get_LRisochron_L,
                                yd=x2calc,hessian=TRUE)
        if (testfit$value < misfit){
            gof <- testfit$value
            bestfit <- testfit
        }
        fit <- list(a=c(testfit$par[3],0.01),
                    b=c(testfit$par[4],0.01),
                    cov.ab=0)
        scatterplot(x2calc,fit=fit)
    }
    covmat <- inverthess(bestfit$hessian)
    a <- sgn*bestfit$par[3]
    sa <- sqrt(covmat[3,3])
    b <- sgn*bestfit$par[4]
    sb <- sqrt(covmat[4,4])
    list(a=c('a'=unname(a),'s[a]'=unname(sa)),
         b=c('b'=unname(b),'s[b]'=unname(sb)),
         cov.ab=unname(covmat[3,4]),
         xyz=x,model=4,n=nrow(x2calc))
}
#' @param anchor control parameters to fix the intercept age or
#'     non-radiogenic composition of the isochron fit. This can be a
#'     scalar or a vector.
#'
#' If \code{anchor[1]=0}: do not anchor the isochron.
#'
#' If \code{anchor[1]=1}: fix the non-radiogenic composition at the
#' values stored in \code{settings('iratio',...)}, OR, if \code{x} has
#' class \code{ThU} and \code{x$format} = \code{3} or \code{4}, fix
#' the intercept at the value stored in \code{x$U8Th2}.
#'
#' If \code{anchor[1]=2}: fix the age at the value stored in \code{anchor[2]}.
#' @noRd
LRisochron.PbPb <- function(x,inverse=TRUE,anchor=0,hide=NULL,omit=NULL,...){
    yd <- data2york(x,inverse=TRUE)
    yd2calc <- clear(yd,hide,omit)
    if (anchor[1]<1){
        out <- LRisochron.default(yd2calc,left=TRUE)
    } else if (anchor[1]==1){
        Pb74 <- iratio('Pb207Pb204')[1]
        Pb64 <- iratio('Pb206Pb204')[1]
        y0 <- Pb74/Pb64
        offset <- 1/Pb64
        yd2calc[,'X'] <- yd2calc[,'X'] - offset
        fit <- anchoredLRisochron(yd2calc,fun=yd2ratios_left,y0=y0)
        stop("TODO")
    } else if (anchor[1]==2 & length(anchor)>1){
        y0 <- age2ratio(tt=anchor[2],ratio="Pb207Pb206")[1]
        out <- anchoredLRisochron(yd2calc,fun=yd2ratios_left,y0=y0)
    } else {
        stop("Invalid anchor.")
    }
    if (inverse){
        out$xyz <- yd
    } else {
        out <- invertfit(out,type='d')
        out$xyz <- normal2inverse(yd,type='d')
    }
    out
}
#' @noRd
LRisochron.ThU <- function(x,inverse=TRUE,anchor=0,hide=NULL,omit=NULL,...){
    if (x$format<3){
        stop("Rightmost isochrons are only available for ThU formats 3 and 4.")
    }
    yd <- data2york(x,inverse=FALSE)
    if (anchor[1]<1){
        out <- LRisochron.default(yd,left=FALSE,hide=hide,omit=omit)
    } else {
        yd2calc <- clear(yd,hide,omit)
        lpi <- 0
        y0i <- gi <- y0 <- b <- NULL # note: y0 = Th2U8i
        if (anchor[1]==1){
            y0i <- y0 <- 1/x$U8Th2
            g <- (yd2calc[,'Y']-y0i)/(yd2calc[,'X']-y0i)
            lsi <- log(stats::sd(g))
            gi <- mean(g)
            out <- anchoredLRisochron(yd2calc, fun=yd2ratios_ThU,
                                      gi=gi, lpi=lpi, lsi=lsi,
                                      y0=y0, gam=NULL)
        } else if (anchor[1]==2 & length(anchor)>1){
            b <- age2ratio(tt=anchor[2],ratio='Th230U238')[1]
            leftmost <- which.min(yd[,'X'])
            x2 <- yd2calc[leftmost,'X']
            y2 <- yd2calc[leftmost,'Y']
            y0i <- y2 - b * x2
            g <- (yd2calc[,'Y']-y0i)/(yd2calc[,'X']-y0i)
            lsi <- log(stats::sd(g))
            out <- anchoredLRisochron(yd2calc, fun=yd2ratios_ThU,
                                      gi=NULL, lpi=lpi, lsi=lsi,
                                      y0i=y0i, y0=NULL, gam=b)
        } else {
            stop("Invalid anchor")
        }
    }
    if (inverse){
        out <- invertfit(out,type='d')
        out$xyz <- normal2inverse(yd,type='d')
    } else {
        out$xyz <- yd
    }
    out
}

anchoredLRisochron <- function(yd2calc,fun,gam=NULL,y0=NULL,...){
    lpi <- 0
    if (is.null(gam) && !is.null(y0)){
        anchor <- 'a'
        g <- (yd2calc[,'Y']-y0)/yd2calc[,'X']
        gi <- mean(g)
        lsi <- log(stats::sd(g))
        init <- c(gi,lpi,lsi)
    } else if (is.null(y0) && !is.null(gam)){
        anchor <- 'b'
        
        init <- c(lpi,lsi,y0i)
    } else {
        stop("Either gam or y0 must be specified.")
    }
    fit <- stats::optim(par=init,fn=get_LRisochron_L,
                        yd=yd2calc,y0=y0,gam=gam,
                        fun=fun,...,hessian=TRUE)
    covmat <- inverthess(fit$hessian)
    np <- length(fit$par)
    if (anchor=='a'){
        sy0 <- 0
        b <- fit$par[1]
        sb <- sqrt(covmat[1,1])
        J <- rbind(c(-y0,1-b),
                   c(1,0))
        E <- J %*% rbind(c(covmat[1,1],0),c(0,0)) %*% t(J)
        a <- y0*(1-b)
        sa <- sqrt(E[1,1])
        cov.ab <- E[1,2]
    } else { # anchor == 'b'
        y0 <- fit$par[np]
        sy0 <- sqrt(covmat[np,np])
        b <- if (identical(fun,yd2ratios_left)) -gam else gam
        sb <- 0
        a <- y0*(1-b)
        sa <- sy0*(1-b)
        cov.ab <- 0
    }
    list(a=c('a'=unname(a),'s[a]'=unname(sa)),
         b=c('b'=unname(b),'s[b]'=unname(sb)),
         cov.ab=unname(cov.ab),
         model=4,n=nrow(yd2calc))
}

init_LRisochron <- function(yd,a=NULL,b=NULL,left=TRUE){
    X <- yd[,'X']
    Y <- yd[,'Y']
    anchor <- 0
    if (is.null(a) && is.null(b)){
        fit <- stats::lm(Y ~ X)
    } else if (is.null(b)){
        anchor <- 'a'
        fit <- stats::lm(I(Y - a) ~ X - 1)
    } else {
        anchor <- 'b'
        fit <- stats::lm(Y ~ 1 + offset(b * X))
    }
    if (summary(fit)$r.squared>0.8){
        a <- fit$coefficients[1]
        b <- fit$coefficients[2]
    } else {
        if (anchor=='a'){
            if (left){
                b <- max((Y-a)/X)
            } else {
                b <- min((Y-a)/X)
            }
        } else if (anchor=='b'){
            if (left){
                a <- max(Y - b*X)
            } else {
                a <- min(Y - b*X)
            }
        } else { # not anchored
            hull_indices <- chull(X,Y)
            if (left){
                leftmost <- which.min(X)
                topmost <- which.max(Y)
                a <- Y[leftmost]
                b <- (Y[topmost]-a)/(X[topmost]-Y[leftmost])
            } else {
                rightmost <- which.max(X)
                a <- min(Y)
                b <- (Y[rightmost]-a)/X[rightmost]
            }
        }
    }
    gi <- ifelse(left,-b,b)
    lpi <- 0
    lsi <- log(stats::sd((Y-a)/X))
    y0i <- a
    if (anchor=='a'){
        out <- c(gi,lpi,lsi)
    } else if (anchor=='b'){
        out <- c(lpi,lsi,y0i)
    } else {
        out <- c(gi,lpi,lsi,y0i)
    }
    out
}

yd2ratios <- function(yd,y0){
    r <- (yd[,'Y']-y0)/yd[,'X']
    sr <- sqrt(errorprop1x2(J1=-r/yd[,'X'],
                            J2=1/yd[,'X'],
                            E11=yd[,'sX']^2,
                            E22=yd[,'sY']^2,
                            E12=yd[,'rXY']*yd[,'sX']*yd[,'sY']))
    cbind(r,sr)
}

get_LRisochron_L <- function(pars,yd,a=NULL,b=NULL){
    np <- length(pars)
    if (is.null(a)){
        y0 <- pars[3]
    }
    if (is.null(b)){
        gam <- pars[np]
    }
    prop <- 1/(exp(pars[1])+1)
    sig <- exp(pars[2])
    mu <- gam
    zs <- yd2ratios(yd=yd,y0=y0)
    z <- zs[,1]
    s <- zs[,2]
    AA  <- prop/sqrt(2*pi*s^2)
    BB <- -0.5*((z-gam)/s)^2
    CC <- (1-prop)/sqrt(2*pi*(sig^2+s^2))
    mu0 <- (mu/sig^2 + z/s^2)/(1/sig^2 + 1/s^2)
    s0 <- 1/sqrt(1/sig^2 + 1/s^2)
    DD <- 1-stats::pnorm((gam-mu0)/s0)
    EE <- 1-stats::pnorm((gam-mu)/sig)
    FF <- -0.5*((z-mu)^2)/(sig^2+s^2)
    fu <- AA*exp(BB) + CC*(DD/EE)*exp(FF)
    fu[fu<.Machine$double.xmin] <- .Machine$double.xmin
    fu[fu>.Machine$double.xmax] <- .Machine$double.xmax
    sum(-log(fu))
}
