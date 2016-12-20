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
#' @param show.numbers boolean flag (\code{TRUE} to show grain
#'     numbers)
#' @param pch plot character (default is a filled circle)
#' @param bg background colour of the plot character
#' @param title add a title to the plot?
#' @param k number of peaks to fit using the finite mixture models of
#'     Galbraith and Green (1993). Setting \code{k='auto'}
#'     automatically selects an optimal number of components based on
#'     the Bayes Information Criterion (BIC). Setting \code{k='min'}
#'     estimates the minimum value using a three parameter model
#'     consisting of a Normal distribution truncated by a discrete
#'     component.
#' @param markers vector of ages of radial marker lines to add to the
#'     plot.
#' @param ... additional arguments to the generic \code{points}
#'     function
#' @references Galbraith, R.F., 1990. The radial plot: graphical
#'     assessment of spread in ages. International Journal of
#'     Radiation Applications and Instrumentation. Part D. Nuclear
#'     Tracks and Radiation Measurements, 17(3), pp.207-214.
#'
#' Galbraith, R.F. and Laslett, G.M., 1993. Statistical models for
#' mixed fission track ages. Nuclear tracks and radiation
#' measurements, 21(4), pp.459-470.
#' @examples
#' data(examples)
#' radialplot(examples$FT1)
#' @rdname radialplot
#' @export
radialplot <- function(x,...){ UseMethod("radialplot",x) }
#' @rdname radialplot
#' @export
radialplot.default <- function(x,from=NA,to=NA,t0=NA,transformation='log',
                               sigdig=2,show.numbers=FALSE,pch=21,
                               bg='white',k=0,markers=NULL,...){
    peaks <- peakfit(x,k=k,sigdig=sigdig)
    markers <- c(markers,peaks$peaks)
    X <- x2zs(x,t0=t0,from=from,to=to,transformation=transformation)
    radial.plot(X,show.numbers=show.numbers,pch=pch,bg=bg,markers=markers,...)
    if (!is.null(peaks$legend))
        graphics::legend('bottomleft',legend=peaks$legend,bty='n')
}
#' @param exterr propagate the external sources of uncertainty into
#'     the mixture model errors?
#' @rdname radialplot
#' @export
radialplot.fissiontracks <- function(x,from=NA,to=NA,t0=NA,
                                     transformation='arcsin',
                                     sigdig=2,show.numbers=FALSE,
                                     pch=21,bg='white',title=TRUE,
                                     markers=NULL,k=0,exterr=TRUE,...){
    peaks <- peakfit(x,k=k,exterr=exterr,sigdig=sigdig)
    markers <- c(markers,peaks$peaks)
    X <- x2zs(x,t0=t0,from=from,to=to,transformation=transformation)
    radial.plot(X,zeta=x$zeta[1],rhoD=x$rhoD[1],
                show.numbers=show.numbers,pch=pch,
                bg=bg,markers=markers,...)
    if (title) title(radial.title(central(x),sigdig=sigdig))
    if (!is.null(peaks$legend))
        graphics::legend('bottomleft',legend=peaks$legend,bty='n')
}
#' @param type scalar indicating whether to plot the
#'     \eqn{^{207}}Pb/\eqn{^{235}}U age (\code{type}=1), the
#'     \eqn{^{206}}Pb/\eqn{^{238}}U age (\code{type}=2), the
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb age (type=3), the
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb-\eqn{^{206}}Pb/\eqn{^{238}}U age
#'     (\code{type}=4), or the (Wetherill) concordia age
#'     (\code{type}=5)
#' @param cutoff.76 the age (in Ma) below which the
#'     \eqn{^{206}}Pb/\eqn{^{238}}U and above which the
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb age is used. This parameter is
#'     only used if \code{type=4}.
#' @param cutoff.disc two element vector with the maximum and minimum
#'     percentage discordance allowed between the
#'     \eqn{^{207}}Pb/\eqn{^{235}}U and \eqn{^{206}}Pb/\eqn{^{238}}U
#'     age (if \eqn{^{206}}Pb/\eqn{^{238}}U < \code{cutoff.76}) or
#'     between the \eqn{^{206}}Pb/\eqn{^{238}}U and
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb age (if
#'     \eqn{^{206}}Pb/\eqn{^{238}}U > \code{cutoff.76}).  Set
#'     \code{cutoff.disc=NA} if you do not want to use this filter.
#' @rdname radialplot
#' @export
radialplot.UPb <- function(x,from=NA,to=NA,t0=NA,
                           transformation='log',type=4,
                           cutoff.76=1100, cutoff.disc=c(-15,5),
                           show.numbers=FALSE, pch=21,bg='white',
                           markers=NULL,k=0,exterr=TRUE,...){
    peaks <- peakfit(x,k=k,exterr=exterr)
    markers <- c(markers,peaks$peaks)
    age2radial(x,from=from,to=to,t0=t0,transformation=transformation,
               type=type,cutoff.76=cutoff.76,cutoff.disc=cutoff.disc,
               show.numbers=show.numbers,pch=pch,bg=bg,markers=markers,
               k=k,...)
    if (!is.null(peaks$legend))
        graphics::legend('bottomleft',legend=peaks$legend,bty='n')
}
#' @rdname radialplot
#' @export
radialplot.ArAr <- function(x,from=NA,to=NA,t0=NA,
                            transformation='log', show.numbers=FALSE,
                            pch=21,bg='white',markers=NULL,k=0,
                            exterr=TRUE,...){
    peaks <- peakfit(x,k=k,exterr=exterr)
    markers <- c(markers,peaks$peaks)
    age2radial(x,from=from,to=to,t0=t0,transformation=transformation,
               show.numbers=show.numbers,pch=pch,bg=bg,markers=markers,
               k=k,...)
    if (!is.null(peaks$legend))
        graphics::legend('bottomleft',legend=peaks$legend,bty='n')
}
#' @rdname radialplot
#' @export
radialplot.UThHe <- function(x,from=NA,to=NA,t0=NA,
                             transformation='log',show.numbers=FALSE,
                             pch=21,bg='white',markers=NULL,k=0,...){
    peaks <- peakfit(x,k=k)
    markers <- c(markers,peaks$peaks)
    age2radial(x,from=from,to=to,t0=t0,transformation=transformation,
               show.numbers=show.numbers,pch=pch,bg=bg,markers=markers,
               k=k,...)
    if (!is.null(peaks$legend))
        graphics::legend('bottomleft',legend=peaks$legend,bty='n')
}

age2radial <- function(x,from=NA,to=NA,t0=NA,transformation='log',
                       type=4,cutoff.76=1100,cutoff.disc=c(-15,5),
                       show.numbers=FALSE,pch=21,bg='white',
                       markers=NULL,k=0,...){
    if (hasClass(x,'UPb')){
        tt <- filter.UPb.ages(x,type,cutoff.76,
                              cutoff.disc,exterr=FALSE)
    } else if (hasClass(x,'ArAr')){
        tt <- ArAr.age(x,exterr=FALSE)
    } else if (hasClass(x,'UThHe')){
        tt <- UThHe.age(x)
    }
    radialplot.default(tt,from=from,to=to,t0=t0,
                       transformation=transformation,
                       show.numbers=show.numbers,pch=pch,bg=bg,
                       markers=markers,...)
}

radial.plot <- function(x,zeta=0,rhoD=0,asprat=3/4,
                        show.numbers=FALSE, pch=21,bg='white',
                        markers=NULL,...){
    exM <- radial.scale(x,zeta,rhoD)
    tticks <- get.radial.tticks(x)
    plot.radial.lines(tticks,l=0.025,x,exM[1],exM[2],
                      zeta,rhoD,label=TRUE)
    if (!is.null(markers)){
        plot.radial.lines(markers,x,exM[1],exM[2],
                          zeta,rhoD,label=FALSE)
    }
    plot.radial.axes(x)
    plot.points(x,show.numbers,pch,bg,...)
}

plot.points <- function(x,show.numbers=FALSE,pch=21,bg='white',...){
    rxy <- data2rxry(x)
    rx <- rxy[,1]
    ry <- rxy[,2]
    if (show.numbers) {
        if('cex' %in% names(list(...)))
            points(rx,ry,pch=pch,bg=bg,...)
        else
            points(rx,ry,pch=pch,bg=bg,cex=3,...)
        text(rx,ry,1:length(rx))
    } else {
        points(rx,ry,pch=pch,bg=bg,...)
    }
}

plot.radial.axes <- function(x){
    xs <- stats::na.omit(x$s)
    graphics::Axis(side=2,at=c(-2,0,2),labels=c(-2,0,2))
    if (identical(x$transformation,'arcsin')){
        plabels <- pretty(c(0,range(1/(2*xs)^2 - 1/2)))
        pticks <- (2*sqrt(plabels+1/2))
    } else {
        plabels <- pretty(c(0,1/xs))
        pticks <- plabels
    }
    graphics::Axis(side=1,at=pticks,labels=plabels)
}

data2rxry <- function(x){
    rx <- 1/x$s
    ry <- (x$z-x$z0)/x$s
    cbind(rx,ry)
}

radial.scale <- function(x,zeta=0,rhoD=0){
    zm <- t2z(x$from,x,zeta,rhoD)
    zM <- t2z(x$to,x,zeta,rhoD)
    padding <- 1.1
    N <- 50
    a <- grDevices::dev.size()[2]/grDevices::dev.size()[1] # aspect ratio
    e <- a/(zM-zm) # ellipticity of the arc
    # get rxM
    theta <- atan(e*(x$z-x$z0))
    rxy <- data2rxry(x)
    xM <- padding * sqrt(max(
          (rxy[,1]^2+rxy[,2]^2)/(cos(theta)^2+(sin(theta)/e)^2), na.rm=TRUE
          ))
    # plot arc
    Z <- seq(zm-x$z0,zM-x$z0,length.out=N)
    rxy <- z2rxy(Z,e,xM)
    graphics::plot(rxy[,1],rxy[,2],type='l',xlim=c(0,xM),axes=FALSE,bty='n',
                   xlab=x$xlab,ylab='standardised estimate')
    c(e,xM)
}

plot.radial.lines <- function(tt,x,e,xM,zeta=0,rhoD=0,l=1,label=FALSE){
    z <- t2z(tt,x,zeta,rhoD)
    rxyb <- z2rxy(z-x$z0,e,xM)
    rxye <- z2rxy(z-x$z0,e,(1-l)*xM)
    for (i in 1:length(tt)){
        graphics::lines(c(rxyb[i,1],rxye[i,1]),
                        c(rxyb[i,2],rxye[i,2]))
        if (label) {
            graphics::text(rxyb[i,1],rxyb[i,2],
                           labels=tt[i],pos=4,xpd=NA)
        }
    }
}

z2rxy <- function(Z,e,xM){
    theta <- atan(e*Z)
    rx <- xM*cos(theta)
    ry <- (xM/e)*sin(theta)
    cbind(rx,ry)
}

t2z <- function(tt,x,zeta,rhoD){
    if (identical(x$transformation,'log')){
        out <- log(tt+x$offset)
    } else if (identical(x$transformation,'arcsin')){
        out <- att(tt,zeta,rhoD)
    } else if (identical(x$transformation,'linear')){
        out <- tt
    }
    out
}

get.radial.tticks <- function(x){
    if (identical(x$transformation,'linear')){
        out <- pretty(c(x$from,x$to))
    } else if (identical(x$transformation,'log')){
        logrange <- log10(c(x$from,x$to)+x$offset)
        out <- grDevices::axisTicks(usr=logrange,log=TRUE)-x$offset
    } else {
        logrange <- log10(c(x$from,x$to))
        out <- grDevices::axisTicks(usr=logrange,log=TRUE)
    }
    nt <- length(out)
    reldiff <- (x$to-out[nt])/(out[nt]-out[nt-1])
    if (reldiff > 0.25) {
        sigdig <- ceiling(1-log10(1-out[nt]/x$to))
        out <- c(out,signif(x$to,sigdig))
    }
    reldiff <- (out[1]-x$from)/(out[2]-out[1])
    if (reldiff > 0.25) {
        sigdig <- ceiling(1-log10(abs(1-x$from/out[1])))
        out <- c(signif(x$from,sigdig),out)
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
x2zs.default <- function(x,t0=NA,from=NA,to=NA,transformation=NA){
    out <- list()
    if (is.na(transformation)) out$transformation <- 'log'
    else out$transformation <- transformation
    if (identical(transformation,'log')){
        out$offset <- get.offset(x[,1],from)
        out$z <- log(x[,1]+out$offset)
        out$s <- x[,2]/(x[,1]+out$offset)
        if (out$offset>0){
            out$xlab <- substitute('t/('~sigma~'+'~a~')',list(a=out$offset))
        } else {
            out$xlab <- expression(t/sigma)
        }
    } else {
        out$z <- x[,1]
        out$s <- x[,2]
        out$xlab <- expression(1/sigma)
    }
    out$z0 <- get.z0(out,t0,from,to)
    # reset limits if necessary
    if (is.na(from)){
        if (identical(transformation,'log'))
            out$from <- exp(min(out$z,na.rm=TRUE))-out$offset
        else
            out$from <- min(out$z,na.rm=TRUE)
    } else {
        out$from <- from
    }
    if (is.na(to)){
        if (identical(transformation,'log'))
            out$to <- exp(max(out$z,na.rm=TRUE))-out$offset
        else if (identical(transformation,'linear'))
            out$to <- max(out$z,na.rm=TRUE)
    } else {
        out$to <- to
    }
    out$zlim <- range(out$z,na.rm=TRUE)
    out
}
x2zs.fissiontracks <- function(x,t0=NA,from=NA,to=NA,transformation=NA){
    out <- list()
    if (is.na(transformation)){
        if (x$format==1) out$transformation <- 'arcsin'
        else out$transformation <- 'log'
    } else {
        out$transformation <- transformation
    }
    if (x$format==1){
        Ns <- x$x[,'Ns']
        Ni <- x$x[,'Ni']
        if (identical(transformation,'linear')){
            tt <- fissiontrack.age(x,exterr=FALSE)
            out$z <- tt[,'t']
            out$s <- tt[,'s[t]']
            out$z0 <- get.z0(out,t0,from,to)
            out$xlab <- expression(1/sigma)
        }
        if (identical(transformation,'log')){
            tt <- fissiontrack.age(x,exterr=FALSE)
            if (any(tt[,'t']<=0)) {
                out$transformation <- 'arcsin'
            } else {
                out$offset <- 0
                out$z <- log(tt[,'t'])
                out$s <- tt[,'s[t]']/tt[,'t']
                out$z0 <- get.z0(out,t0,from,to)
                out$xlab <- expression(t/sigma)
            }
        }
        if (identical(out$transformation,'arcsin')){
            out$z <- atan(sqrt((Ns+3/8)/(Ni+3/8)))
            out$s <- sqrt(1/(Ns+Ni+1/2))/2
            if (is.na(t0))
                if (is.na(from) | is.na(to)) {
                    out$z0 <- atan(sqrt(sum(Ns,na.rm=TRUE)/
                                        sum(Ni,na.rm=TRUE)))
                } else {
                    zmin <- att(from,x$zeta[1],x$rhoD[1])
                    zmax <- att(to,x$zeta[1],x$rhoD[1])
                    out$z0 <- mean(c(zmin,zmax),na.rm=TRUE)
                }
            else {
                out$z0 <- att(t0,x$zeta[1],x$rhoD[1])
            }
            out$xlab <- 'Ns+Ni'
        }
        # reset limits if necessary
        out$zlim <- range(out$z,na.rm=TRUE)
        if (is.na(from)){
            if (identical(transformation,'log'))
                out$from <- exp(min(out$z,na.rm=TRUE))
            else if (identical(transformation,'arcsin'))
                out$from <- iatt(min(out$z,na.rm=TRUE),
                                 x$zeta[1],x$rhoD[1])
            else if (identical(transformation,'linear'))
                out$from <- min(out$z,na.rm=TRUE)
        } else {
            out$from <- from
        }
        if (is.na(to)){
            if (identical(transformation,'log'))
                out$to <- exp(max(out$z,na.rm=TRUE))
            else if (identical(transformation,'arcsin'))
                out$to <- iatt(max(out$z,na.rm=TRUE),
                               x$zeta[1],x$rhoD[1])
            else if (identical(transformation,'linear'))
                out$to <- max(out$z,na.rm=TRUE)
        } else {
            out$to <- to
        }
    } else {
        tt <- fissiontrack.age(x,exterr=FALSE)
        if (identical(transformation,'arcsin')) transformation <- 'log'
        out <- x2zs.default(tt,t0=t0,from=from,to=to,transformation=transformation)
    }
    out
}

# only for log and lin transformation
# x is a list containing the items z, s, transformation
# and (if transformation=='log) offset
get.z0 <- function(x,t0=NA,from=NA,to=NA){
    if (is.na(t0)){
        if (is.na(from) | is.na(to)) {
            z0 <- mean(x$z,na.rm=TRUE)
        } else if (identical(x$transformation,'log')) {
            z0 <- mean(log(c(from,to)+x$offset),na.rm=TRUE)
        } else if (identical(x$transformation,'linear')) {
            z0 <- mean(c(from,to),na.rm=TRUE)
        } else {
            stop('illegal input')
        }
    } else if (identical(x$transformation,'log')){
        z0 <- log(t0+x$offset)
    } else if (identical(x$transformation,'linear')){
        z0 <- t0
    } else {
        stop('illegal input')
    }
    z0
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

radial.title <- function(fit,sigdig=2){
    rounded.age <- roundit(fit$age[1],fit$age[2],sigdig=sigdig)
    line1 <- substitute('central age ='~a%+-%b~'(1'~sigma~')',
                        list(a=rounded.age$x, b=rounded.age$err))
    line2 <- substitute('dispersion ='~a~'%, p('~chi^2*')='~b,
                        list(a=signif(100*fit$disp,2),
                             b=signif(fit$p.value,2)))
    graphics::mtext(line1,line=1)
    graphics::mtext(line2,line=0)
}

get.offset <- function(x,from=NA){
    m <- min(c(x,from),na.rm=TRUE)
    if (m>0){
        offset <- 0;
    } else if (m==0){
        offset <- 1;
    } else {
        offset = 10^(floor(log10(-m))+1);
    }
    offset
}
