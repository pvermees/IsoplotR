#' Visualise heteroscedastic data on a radial plot
#'
#' Implementation of a graphical device developed by Rex Galbraith to
#' display several estimates of the same quantity that have different
#' standard errors.
#'
#' @param x Either an nx2 matix of (transformed) values z and their
#'     standard errors s
#'
#' OR
#'
#' and object of class \code{fissiontracks}, \code{UThHe},
#' \code{ArAr}, or \code{UPb}
#' @param from minimum age limit of the radial scale
#' @param to maximum age limit of the radial scale
#' @param t0 central value
#' @param transformation one of either \code{log}, \code{linear} or
#'     (if \code{x} has class \code{fissiontracks})
#' @param sigdig the number of significant digits of the numerical
#'     values reported in the title of the graphical output.
#' @param show.numbers boolean flag (TRUE to show grain numbers)
#' @param pch plot character (default is a filled circle)
#' @param bg background colour of the plot character
#' @param ... additional arguments to the generic \code{points} function
#' @examples
#' data(examples)
#' radialplot(examples$FT1)
#' @rdname radialplot
#' @export
radialplot <- function(x,...){ UseMethod("radialplot",x) }
#' @rdname radialplot
#' @export
radialplot.default <- function(x,from=NA,to=NA,t0=NA,
                               transformation='log',sigdig=2,
                               show.numbers=FALSE,pch=21,
                               bg='white',...){
    X <- x2zs(x,transformation=transformation,t0=t0)
    radial.plot(X,from=from,to=to,
                transformation=transformation,
                show.numbers=show.numbers,
                pch=pch,bg=bg,...)
}
#' @rdname radialplot
#' @export
radialplot.fissiontracks <- function(x,from=NA,to=NA,t0=NA,
                                     transformation='arcsin',
                                     sigdig=2,show.numbers=FALSE,
                                     pch=21,bg='white',...){
    if (x$format>1 & transformation=='arcsin') transformation <- 'log'
    X <- x2zs(x,transformation=transformation,t0=t0)
    radial.plot(X,from=from,to=to,zeta=x$zeta[1],
                rhoD=x$rhoD[1],transformation=transformation,
                show.numbers=show.numbers,pch=pch,bg=bg,...)
}
#' @param type scalar indicating whether to plot the
#'     \eqn{^{207}}Pb/\eqn{^{235}}U age (type=1), the
#'     \eqn{^{206}}Pb/\eqn{^{238}}U age (type=2), the
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb age (type=3), the
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb-\eqn{^{206}}Pb/\eqn{^{238}}U age
#'     (type=4), or the (Wetherill) concordia age (type=5)
#' @param cutoff.76 the age (in Ma) below which the
#'     \eqn{^{206}}Pb/\eqn{^{238}}U and above which the
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb age is used. This parameter is
#'     only used if \code{type=4}.
#' @param cutoff.disc two element vector with the maximum and minimum
#'     percentage discordance allowed between the
#'     \eqn{^{207}}Pb/\eqn{^{235}}U and \eqn{^{206}}Pb/\eqn{^{238}}U
#'     age (if \eqn{^{206}}Pb/\eqn{^{238}}U < cutoff.76) or between
#'     the \eqn{^{206}}Pb/\eqn{^{238}}U and
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb age (if
#'     \eqn{^{206}}Pb/\eqn{^{238}}U > cutoff.76).  Set
#'     \code{cutoff.disc=NA} if you do not want to use this filter.
#' @rdname radialplot
#' @export
radialplot.UPb <- function(x,from=NA,to=NA,t0=NA,
                           transformation='log',
                           type=4,cutoff.76=1100,
                           cutoff.disc=c(-15,5),
                           sigdig=2,show.numbers=FALSE,
                           pch=21,bg='white',...){
    tt <- filter.UPb.ages(x,type,cutoff.76,
                          cutoff.disc,exterr=FALSE)
    radialplot.default(tt,from=from,to=to,t0=t0,
                       transformation=transformation,
                       show.numbers=show.numbers,
                       pch=pch,bg=bg,...)
}
#' @rdname radialplot
#' @export
radialplot.ArAr <- function(x,from=NA,to=NA,t0=NA,
                            transformation='log',
                            sigdig=2,show.numbers=FALSE,
                            pch=21,bg='white',...){
    tt <- ArAr.age(x,exterr=FALSE)
    radialplot.default(tt,from=from,to=to,t0=t0,
                       transformation=transformation,
                       show.numbers=show.numbers,
                       pch=pch,bg=bg,...)
}
#' @rdname radialplot
#' @export
radialplot.UThHe <- function(x,from=NA,to=NA,t0=NA,
                             transformation='log',
                             sigdig=2,show.numbers=FALSE,
                             pch=21,bg='white',...){
    tt <- UThHe.age(x)
    radialplot.default(tt,from=from,to=to,t0=t0,
                       transformation=transformation,
                       show.numbers=show.numbers,
                       pch=pch,bg=bg,...)
}

radial.plot <- function(x,from=NA,to=NA,zeta=0,rhoD=0,
                        transformation='arcsin',asprat=3/4,
                        show.numbers=FALSE,
                        pch=21,bg='white',...){
    # 1. get z and t limits of the radial scale
    zlim <- range(x$z) # initialise using linear transformation
    tlim <- c(from,to)
    if (is.na(from)){
        if (identical(transformation,'log'))
            from <- exp(min(x$z))
        else if (identical(transformation,'arcsin'))
            from <- iatt(min(x$z),zeta,rhoD)
        else if (identical(transformation,'linear'))
            from <- min(x$z)
    }
    if (is.na(to)){
        if (identical(transformation,'log'))
            to <- exp(zlim[2])
        else if (identical(transformation,'arcsin'))
            to <- iatt(zlim[2],zeta,rhoD)
        else if (identical(transformation,'linear'))
            to <- zlim[2]
    }
    # 2. get time ticks and z limits
    tticks <- get.radial.tticks(from,to,transformation)
    if (identical(transformation,'log'))
        zticks <- log(tticks)
    else if (identical(transformation,'arcsin'))
        zticks <- att(tticks,zeta,rhoD)
    else if (identical(transformation,'linear'))
        zticks <- tticks
    zlim <- range(zticks)
    # 3. get axis labels
    if (identical(transformation,'log'))
        xlab <- expression(t/sigma)
    else if (identical(transformation,'arcsin'))
        xlab <- 'Ns+Ni'
    else
        xlab <- expression(1/sigma)
    ylab <- 'standardised estimate'
    # 4. draw arc
    fz <- stats::optim(0.5,get.fz,method='BFGS',
                       z0=x$z0,zlim=zlim,asprat=asprat)$par
    fxy <- get.fxy(x,fz,asprat)
    zscale <- seq(zlim[1],zlim[2],length.out=50)
    rx <- cos(fz*(zscale-x$z0))
    ry <- sin(fz*(zscale-x$z0))
    graphics::plot(rx,ry,type='l',xlim=c(0,1),bty='n',
                   axes=FALSE,xlab=xlab,ylab=ylab)
    # 5. draw the ticks
    for (i in 1:length(tticks)){
        rxb <- cos(fz*(zticks[i]-x$z0))
        ryb <- sin(fz*(zticks[i]-x$z0))
        rxe <- 0.98*rxb
        rye <- 0.98*ryb
        graphics::lines(c(rxb,rxe),c(ryb,rye))
        graphics::text(rxb,ryb,labels=tticks[i],
                       pos=4,xpd=NA)
    }
    # 6. plot points
    rx <- fxy/x$s
    ry <- fxy*fz*(x$z-x$z0)/x$s
    if (show.numbers) {
        if('cex' %in% names(list(...)))
            points(rx,ry,pch=pch,bg=bg,...)
        else
            points(rx,ry,pch=pch,bg=bg,cex=3,...)
        text(rx,ry,1:length(rx))
    } else {
        points(rx,ry,pch=pch,bg=bg,...)
    }
    # 7. x and y axes
    graphics::Axis(side=2,at=fxy*fz*c(-2,0,2),labels=c(-2,0,2))
    if (identical(transformation,'arcsin')){
        plabels <- pretty(c(0,range(1/(2*x$s)^2 - 1/2)))
        pticks <- fxy*(2*sqrt(plabels+1/2))
    } else {
        plabels <- pretty(c(0,1/x$s))
        pticks <- fxy*plabels
    }
    graphics::Axis(side=1,at=pticks,labels=plabels)
}

get.radial.tticks <- function(from,to,transformation){
    if (identical(transformation,'linear'))
        out <- pretty(c(from,to))
    else
        out <- grDevices::axisTicks(usr=log10(c(from,to)),log=TRUE)
    nt <- length(out)
    reldiff <- (to-out[nt])/(out[nt]-out[nt-1])
    if (reldiff > 0.25) {
        sigdig <- ceiling(1-log10(1-out[nt]/to))
        out <- c(out,signif(to,sigdig))
    }
    reldiff <- (out[1]-from)/(out[2]-out[1])
    if (reldiff > 0.25) {
        sigdig <- ceiling(1-log10(1-from/out[1]))
        out <- c(signif(from,sigdig),out)
    }
    out
}

get.fxy <- function(x,fz,asprat,buffer=0.9){
    xx <- 1/x$s
    yy <- fz*(x$z-x$z0)/x$s
    rd <- sqrt(xx^2+yy^2) # radial distance
    buffer/max(rd)
}

get.fz <- function(fz,z0,zlim,asprat){
    (sin(atan(fz*(zlim[2]-z0))) + sin(atan(fz*(z0-zlim[1]))) - asprat)^2
}

x2zs <- function(x,...){ UseMethod("x2zs",x) }
x2zs.default <- function(x,transformation='log',t0=NA,...){
    out <- list()
    if (identical(transformation,'log')){
        out$z <- log(x[,1])
        out$s <- x[,2]/x[,1]
        out$xlab <- expression(t/sigma)
        if (is.na(t0)) out$z0 <- mean(out$z)
        else out$z0 <- log(t0)
    } else {
        out$z <- x[,1]
        out$s <- x[,2]
        out$xlab <- expression(1/sigma)
        if (is.na(t0)) out$z0 <- mean(out$z)
        else out$z0 <- t0
    }
    out
}
x2zs.fissiontracks <- function(x,transformation='arcsin',t0=NA,...){
    out <- list()
    if (x$format==1){
        Ns <- x$x[,'Ns']
        Ni <- x$x[,'Ni']
        if (identical(transformation,'linear')){
            tt <- fissiontrack.age(x,exterr=FALSE)
            out$z <- tt[,'t']
            out$s <- tt[,'s[t]']
            if (is.na(t0)) out$z0 <- mean(out$z)
            else out$z0 <- mean(t0)
        }
        if (identical(transformation,'log')){
            tt <- fissiontrack.age(x,exterr=FALSE)
            if (any(tt[,'t']<=0)) {
                transformation <- 'arcsin'
            } else {
                out$z <- log(tt[,'t'])
                out$s <- tt[,'s[t]']/tt[,'t']
                if (is.na(t0)) out$z0 <- mean(out$z)
                else out$z0 <- log(t0)
            }
        }
        if (identical(transformation,'arcsin')){
            out$z <- atan(sqrt((Ns+3/8)/(Ni+3/8)))
            out$s <- sqrt(1/(Ns+Ni+1/2))/2
            if (is.na(t0))
                out$z0 <- atan(sqrt(sum(Ns)/sum(Ni)))
            else
                out$z0 <- att(t0,x$zeta[1],x$rhoD[1])
        }
    } else {
        tt <- fissiontrack.age(x,exterr=FALSE)
        out <- x2zs.default(tt,transformation=transformation,t0=t0)
    }
    out
}

# arctan transformation
att <- function(tt,zeta,rhoD){
    L8 <- lambda('U238')[1]
    atan(sqrt((exp(L8*tt)-1)/(L8*(zeta/1e6)*rhoD/2)))
}
# inverse arctan transformation
iatt <- function(z,zeta,rhoD){
    L8 <- lambda('U238')[1]
    log(1+L8*(zeta/2e6)*rhoD*tan(z)^2)/L8
}
