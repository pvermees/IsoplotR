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
#' and object of class \code{fissiontracks}
#'
#' @rdname radialplot
#' @export
radialplot <- function(x,...){ UseMethod("radialplot",x) }
#' @rdname radialplot
#' @export
radialplot.default <- function(x,from=NA,to=NA,z0=NA,transformation='log',...){
    if (identical(transformation,'log')){
        z <- log(x[,1])
        s <- x[,2]/x[,1]
        xlab <- expression(1/sigma)
    } else {
        z <- x[,1]
        s <- x[,2]
        xlab <- expression(1/sigma)
    }
    z0 <- mean(z)
    rs <- radial.scale(z,s,z0,from=from,to=to,
                       transformation=transformation,...)
    radial.plot(z,s,z0,rs,...)
}
#' @rdname radialplot
#' @export
radialplot.fissiontracks <- function(x,from=NA,to=NA,t0=NA,
                                     transformation='arcsin',...){
    Ns <- x$x[,'Ns']
    Ni <- x$x[,'Ni']
    if (identical(transformation,'arcsin')){
        z <- atan(sqrt((Ns+3/8)/(Ni+3/8)))
        s <- sqrt(1/(Ns+Ni+1/2))/2
        if (is.na(t0))
            z0 <- atan(sqrt(sum(Ns)/sum(Ni)))
        else
            z0 <- att(t0,x$zeta[1],x$rhoD[1])
    } else if (identical(transformation,'log')){
        tt <- fissiontrack.age(x,external=FALSE)
        z <- log(tt[,'t'])
        s <- tt[,'st']/tt[,'t']
        z0 <- mean(z)
    } else if (identical(transformation,'linear')){
        tt <- fissiontrack.age(x,external=FALSE)
        z <- tt[,'t']
        s <- tt[,'st']
        z0 <- mean(z)
    }
    rs <- radial.scale(z,s,z0,from=from,to=to,zeta=x$zeta[1],
                       rhoD=x$rhoD[1],transformation=transformation,...)
    radial.plot(z,s,z0,rs,...)
}

radial.plot <- function(z,s,z0,rs,...){
    x <- list(rx=1/s,ry=(z-z0)/s)    
    xlim <- c(0,max(rs$rx))
    ylim <- c(min(rs$ry),max(rs$ry))
    plot(xlim,ylim,bty='n',axes=FALSE,type='n',
         xlab=rs$xlab,ylab='standardised estimate',
         asp=1,...)
    Axis(side=2,at=c(-2,0,2))
    Axis(side=1,at=rs$xticks,labels=rs$xlabels)
    points(x$rx,x$ry)
    lines(rs$rx,rs$ry)
}

radial.scale <- function(z,s,z0,from=NA,to=NA,zeta=0,rhoD=0,
                         transformation='arcsin',f=1.1){
    zlim <- range(z)
    out <- list()
    nn <- 100
    if (is.na(from)) zlim[1] <- min(z)
    if (is.na(to)) zlim[2] <- max(z)
    if (identical(transformation,'log')){
        if (!is.na(from)) zlim[1] <- log(from)
        if (!is.na(to)) zlim[2] <- max(z)
        tlim <- exp(zlim)
        out$tticks <- pretty(tlim)
        zticks <- log(out$tticks)
        out$xlab <- expression(t/sigma)
        out$xticks <- pretty(range(1/s))
        out$xlabels <- out$xticks
    } else if (identical(transformation,'lin')){
        if (!is.na(from)) zlim[1] <- from
        if (!is.na(to)) zlim[2] <- to
        tlim <- zlim
        out$tticks <- pretty(tlim)
        zticks <- out$tticks
        out$xlab <- expression(1/sigma)
        out$xticks <- pretty(range(1/s))
        out$xlabels <- out$xticks
    } else if (identical(transformation,'arcsin')){
        if (!is.na(from)) zlim[1] <- att(from,zeta,rhoD)
        if (!is.na(to)) zlim[2] <- att(to,zeta,rhoD)
        tlim <- iatt(zlim,zeta,rhoD)
        out$tticks <- pretty(tlim)
        zticks <- att(out$tticks,zeta,rhoD)
        out$xlab <- 'Ns+Ni'
        xrange <- range(1/s)
        NsNiRange <- (xrange/2)^2 - 1/2
        out$xlabels <- pretty(NsNiRange)
        out$xticks <- 2*sqrt(out$xlabels+0.5)
    }
    zrs <- seq(zlim[1],zlim[2],length.out=nn)
    r <- f*max(1/s)
    out$rx <- r*cos(atan(zrs-z0))
    out$ry <- r*sin(atan(zrs-z0))
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
