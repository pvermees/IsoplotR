#' Visualise heteroscedastic data on a radial plot
#'
#' Implementation of a graphical device developed by Rex Galbraith to
#' display several estimates of the same quantity that have different
#' standard errors.
#'
#' @details
#'
#' The radial plot (Galbraith, 1988, 1990) is a graphical device that
#' was specifically designed to display heteroscedastic data, and is
#' constructed as follows.  Consider a set of dates
#' \eqn{\{t_1,...,t_i,...,t_n\}} and uncertainties
#' \eqn{\{s[t_1],...,s[t_i],...,s[t_n]\}}. Define \eqn{z_i = z[t_i]}
#' to be a transformation of \eqn{t_i} (e.g., \eqn{z_i = log[t_i]}),
#' and let \eqn{s[z_i]} be its propagated analytical uncertainty
#' (i.e., \eqn{s[z_i] = s[t_i]/t_i} in the case of a logarithmic
#' transformation). Create a scatterplot of \eqn{(x_i,y_i)} values,
#' where \eqn{x_i = 1/s[z_i]} and \eqn{y_i = (z_i-z_\circ)/s[z_i]},
#' where \eqn{z_\circ} is some reference value such as the mean. The
#' slope of a line connecting the origin of this scatterplot with any
#' of the \eqn{(x_i,y_i)}s is proportional to \eqn{z_i} and, hence,
#' the date \eqn{t_i}.  These dates can be more easily visualised by
#' drawing a radial scale at some convenient distance from the origin
#' and annotating it with labelled ticks at the appropriate
#' angles. While the angular position of each data point represents
#' the date, its horizontal distance from the origin is proportional
#' to the precision. Imprecise measurements plot on the left hand side
#' of the radial plot, whereas precise age determinations are found
#' further towards the right. Thus, radial plots allow the observer to
#' assess both the magnitude and the precision of quantitative data in
#' one glance.
#'
#' @param x Either an \code{[n x 2]} matix of (transformed) values z
#'     and their standard errors s
#'
#' OR
#'
#' and object of class \code{fissiontracks}, \code{UThHe},
#' \code{ArAr}, \code{ReOs}, \code{SmNd}, \code{RbSr}, \code{LuHf},
#' \code{ThU}, \code{PbPb} or \code{UPb}
#' @param from minimum age limit of the radial scale
#' @param to maximum age limit of the radial scale
#' @param t0 central value
#' @param transformation one of either \code{log}, \code{linear} or
#'     (if \code{x} has class \code{fissiontracks}), \code{arcsin}.
#' @param sigdig the number of significant digits of the numerical
#'     values reported in the title of the graphical output.
#' @param show.numbers boolean flag (\code{TRUE} to show grain
#'     numbers)
#' @param pch plot character (default is a filled circle)
#' @param levels a vector with additional values to be displayed as
#'     different background colours of the plot symbols.
#' @param clabel label of the colour legend
#' @param bg a vector of two background colours for the plot symbols.
#'     If \code{levels=NA}, then only the first colour is used. If
#'     \code{levels} is a vector of numbers, then \code{bg} is used to
#'     construct a colour ramp.
#' @param title add a title to the plot?
#' @param k number of peaks to fit using the finite mixture models of
#'     Galbraith and Laslett (1993). Setting \code{k='auto'}
#'     automatically selects an optimal number of components based on
#'     the Bayes Information Criterion (BIC). Setting \code{k='min'}
#'     estimates the minimum value using a three parameter model
#'     consisting of a Normal distribution truncated by a discrete
#'     component.
#' @param markers vector of ages of radial marker lines to add to the
#'     plot.
#' @param alpha cutoff value for confidence intervals
#' @param ... additional arguments to the generic \code{points}
#'     function
#' @seealso \code{\link{peakfit}}, \code{\link{central}}
#' @references
#' Galbraith, R.F., 1988. Graphical display of estimates having
#' differing standard errors. Technometrics, 30(3), pp.271-281.
#'
#' Galbraith, R.F., 1990. The radial plot: graphical assessment of
#' spread in ages. International Journal of Radiation Applications and
#' Instrumentation. Part D. Nuclear Tracks and Radiation Measurements,
#' 17(3), pp.207-214.
#'
#' Galbraith, R.F. and Laslett, G.M., 1993. Statistical models for
#' mixed fission track ages. Nuclear Tracks and Radiation
#' Measurements, 21(4), pp.459-470.
#' @examples
#' data(examples)
#' radialplot(examples$FT1)
#'
#' dev.new()
#' radialplot(examples$LudwigMixture,k='min')
#'
#' @rdname radialplot
#' @export
radialplot <- function(x,...){ UseMethod("radialplot",x) }
#' @rdname radialplot
#' @export
radialplot.default <- function(x,from=NA,to=NA,t0=NA,
                               transformation='log',sigdig=2,
                               show.numbers=FALSE,pch=21,levels=NA,
                               clabel="",bg=c("white","red"),
                               title=TRUE,k=0,markers=NULL,
                               alpha=0.05,...){
    peaks <- peakfit(x,k=k,sigdig=sigdig)
    markers <- c(markers,peaks$peaks['t',])
    X <- x2zs(x,t0=t0,from=from,to=to,transformation=transformation)
    radial.plot(X,show.numbers=show.numbers,pch=pch,levels=levels,
                clabel=clabel,bg=bg,markers=markers,...)
    if (title) title(radial.title(central(x,alpha=alpha),sigdig=sigdig))
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
                                     pch=21,levels=NA,clabel="",
                                     bg=c("white","red"),title=TRUE,
                                     markers=NULL,k=0,exterr=TRUE,
                                     alpha=0.05,...){
    peaks <- peakfit(x,k=k,exterr=exterr,sigdig=sigdig)
    markers <- c(markers,peaks$peaks['t',])
    X <- x2zs(x,t0=t0,from=from,to=to,transformation=transformation)
    radial.plot(X,zeta=x$zeta[1],rhoD=x$rhoD[1],
                show.numbers=show.numbers,pch=pch,levels=levels,
                clabel=clabel,bg=bg,markers=markers,...)
    if (title) title(radial.title(central(x,alpha=alpha),sigdig=sigdig))
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
#' @param common.Pb apply a common lead correction using one of three
#'     methods:
#'
#' \code{1}: use the isochron intercept as the initial Pb-composition
#'
#' \code{2}: use the Stacey-Kramer two-stage model to infer the initial
#' Pb-composition
#'
#' \code{3}: use the Pb-composition stored in
#' \code{settings('iratio','Pb206Pb204')} and
#' \code{settings('iratio','Pb207Pb204')}
#'
#' @rdname radialplot
#' @export
radialplot.UPb <- function(x,from=NA,to=NA,t0=NA,
                           transformation='log',type=4,
                           cutoff.76=1100,cutoff.disc=c(-15,5),
                           show.numbers=FALSE,pch=21,levels=NA,
                           clabel="",bg=c("white","red"),
                           markers=NULL,k=0,exterr=TRUE,common.Pb=0,
                           alpha=0.05,...){
    if (common.Pb %in% c(1,2,3))
        X <- common.Pb.correction(x,option=common.Pb)
    else
        X <- x
    radialplot_helper(X,from=from,to=to,t0=t0,
                      transformation=transformation,type=type,
                      cutoff.76=cutoff.76,cutoff.disc=cutoff.disc,
                      show.numbers=show.numbers,pch=pch,levels=levels,
                      clabel=clabel,bg=bg,markers=markers,k=k,
                      exterr=exterr,alpha=alpha,...)
}
#' @param i2i `isochron to intercept': calculates the initial (aka
#'     `inherited', `excess', or `common')
#'     \eqn{^{40}}Ar/\eqn{^{36}}Ar, \eqn{^{207}}Pb/\eqn{^{204}}Pb,
#'     \eqn{^{87}}Sr/\eqn{^{86}}Sr, \eqn{^{143}}Nd/\eqn{^{144}}Nd,
#'     \eqn{^{187}}Os/\eqn{^{188}}Os or \eqn{^{176}}Hf/\eqn{^{177}}Hf
#'     ratio from an isochron fit. Setting \code{i2i} to \code{FALSE}
#'     uses the default values stored in
#'     \code{settings('iratio',...)}. When applied to data of class
#'     \code{ThU}, setting \code{i2i} to \code{TRUE} applies a
#'     detrital Th-correction.
#' @rdname radialplot
#' @export
radialplot.PbPb <- function(x,from=NA,to=NA,t0=NA,
                            transformation='log',show.numbers=FALSE,
                            pch=21,levels=NA,clabel="",
                            bg=c("white","red"),markers=NULL,k=0,
                            exterr=TRUE,common.Pb=1,alpha=0.05,...){
    if (common.Pb %in% c(1,2,3))
        X <- common.Pb.correction(x,option=common.Pb)
    else
        X <- x
    radialplot_helper(X,from=from,to=to,t0=t0,
                      transformation=transformation,
                      show.numbers=show.numbers,pch=pch,levels=levels,
                      clabel=clabel,bg=bg,markers=markers,k=k,
                      exterr=exterr,alpha=alpha,...)
}
#' @rdname radialplot
#' @export
radialplot.ArAr <- function(x,from=NA,to=NA,t0=NA,
                            transformation='log',show.numbers=FALSE,
                            pch=21,levels=NA,clabel="",
                            bg=c("white","red"),markers=NULL,k=0,
                            exterr=TRUE,i2i=FALSE,alpha=0.05,...){
    radialplot_helper(x,from=from,to=to,t0=t0,
                      transformation=transformation,
                      show.numbers=show.numbers,pch=pch,levels=levels,
                      clabel=clabel,bg=bg,markers=markers,k=k,
                      exterr=exterr,i2i=i2i,alpha=alpha,...)
}
#' @rdname radialplot
#' @export
radialplot.UThHe <- function(x,from=NA,to=NA,t0=NA,
                             transformation='log',show.numbers=FALSE,
                             pch=21,levels=NA,clabel="",
                             bg=c("white","red"),markers=NULL,k=0,
                             alpha=0.05,...){
    radialplot_helper(x,from=from,to=to,t0=t0,
                      transformation=transformation,
                      show.numbers=show.numbers,pch=pch,levels=levels,
                      clabel=clabel,bg=bg,markers=markers,k=k,
                      exterr=FALSE,alpha=alpha,...)
}
#' @rdname radialplot
#' @export
radialplot.ReOs <- function(x,from=NA,to=NA,t0=NA,
                            transformation='log',show.numbers=FALSE,
                            pch=21,levels=NA,clabel="",
                            bg=c("white","red"),markers=NULL,k=0,
                            exterr=TRUE,i2i=TRUE,alpha=0.05,...){
    radialplot_helper(x,from=from,to=to,t0=t0,
                      transformation=transformation,
                      show.numbers=show.numbers,pch=pch,levels=levels,
                      clabel=clabel,bg=bg,markers=markers,k=k,
                      exterr=exterr,i2i=i2i,alpha=alpha,...)
}
#' @rdname radialplot
#' @export
radialplot.SmNd <- function(x,from=NA,to=NA,t0=NA,
                            transformation='log',show.numbers=FALSE,
                            pch=21,levels=NA,clabel="",
                            bg=c("white","red"),markers=NULL,k=0,
                            exterr=TRUE,i2i=TRUE,alpha=0.05,...){
    radialplot_helper(x,from=from,to=to,t0=t0,
                      transformation=transformation,
                      show.numbers=show.numbers,pch=pch,levels=levels,
                      clabel=clabel,bg=bg,markers=markers,k=k,
                      exterr=exterr,i2i=i2i,alpha=alpha,...)
}
#' @rdname radialplot
#' @export
radialplot.RbSr <- function(x,from=NA,to=NA,t0=NA,
                            transformation='log',show.numbers=FALSE,
                            pch=21,levels=NA,clabel="",
                            bg=c("white","red"),markers=NULL,k=0,
                            exterr=TRUE,i2i=TRUE,alpha=0.05,...){
    radialplot_helper(x,from=from,to=to,t0=t0,
                      transformation=transformation,
                      show.numbers=show.numbers,pch=pch,levels=levels,
                      clabel=clabel,bg=bg,markers=markers,k=k,
                      exterr=exterr,i2i=i2i,alpha=alpha,...)
}
#' @rdname radialplot
#' @export
radialplot.LuHf <- function(x,from=NA,to=NA,t0=NA,
                            transformation='log',show.numbers=FALSE,
                            pch=21,levels=NA,clabel="",
                            bg=c("white","red"),markers=NULL,k=0,
                            exterr=TRUE,i2i=TRUE,alpha=0.05,...){
    radialplot_helper(x,from=from,to=to,t0=t0,
                      transformation=transformation,
                      show.numbers=show.numbers,pch=pch,levels=levels,
                      clabel=clabel,bg=bg,markers=markers,k=k,
                      exterr=exterr,i2i=i2i,alpha=alpha,...)
}
#' @rdname radialplot
#' @export
radialplot.ThU <- function(x,from=NA,to=NA,t0=NA,
                           transformation='log',show.numbers=FALSE,
                           pch=21,levels=NA,clabel="",
                           bg=c("white","red"),markers=NULL,k=0,
                           i2i=TRUE,alpha=0.05,...){
    radialplot_helper(x,from=from,to=to,t0=t0,
                      transformation=transformation,
                      show.numbers=show.numbers,pch=pch,levels=levels,
                      clabel=clabel,bg=bg,markers=markers,k=k,
                      exterr=FALSE,i2i=i2i,alpha=alpha,...)
}
radialplot_helper <- function(x,from=NA,to=NA,t0=NA,
                              transformation='log',type=4,
                              cutoff.76=1100,cutoff.disc=c(-15,5),
                              show.numbers=FALSE,pch=21,levels=NA,
                              clabel="",bg=c("white","red"),
                              markers=NULL,k=0,exterr=TRUE,i2i=FALSE,
                              alpha=0.05,...){
    peaks <- peakfit(x,k=k,exterr=exterr,i2i=i2i,type=type,
                     cutoff.76=cutoff.76,cutoff.disc=cutoff.disc)
    markers <- c(markers,peaks$peaks['t',])
    age2radial(x,from=from,to=to,t0=t0,transformation=transformation,
               type=type,cutoff.76=cutoff.76,cutoff.disc=cutoff.disc,
               show.numbers=show.numbers,pch=pch,levels=levels,
               clabel=clabel,bg=bg,markers=markers,i2i=i2i,
               alpha=alpha,...)
    if (!is.null(peaks$legend))
        graphics::legend('bottomleft',legend=peaks$legend,bty='n')
}

age2radial <- function(x,from=NA,to=NA,t0=NA,transformation='log',
                       type=4,cutoff.76=1100,cutoff.disc=c(-15,5),
                       show.numbers=FALSE,pch=21,levels=NA,clabel="",
                       bg=c("white","red"),markers=NULL,k=0,
                       i2i=FALSE,alpha=0.05,...){
    if (hasClass(x,'UPb')){
        tt <- filter.UPb.ages(x,type=type,cutoff.76=cutoff.76,
                              cutoff.disc=cutoff.disc,exterr=FALSE)
    } else if (hasClass(x,'PbPb')){
        tt <- PbPb.age(x,exterr=FALSE,i2i=i2i)
    } else if (hasClass(x,'ArAr')){
        tt <- ArAr.age(x,exterr=FALSE,i2i=i2i)
    } else if (hasClass(x,'UThHe')){
        tt <- UThHe.age(x)
    } else if (hasClass(x,'ReOs')){
        tt <- ReOs.age(x,exterr=FALSE,i2i=i2i)
    } else if (hasClass(x,'SmNd')){
        tt <- SmNd.age(x,exterr=FALSE,i2i=i2i)
    } else if (hasClass(x,'RbSr')){
        tt <- RbSr.age(x,exterr=FALSE,i2i=i2i)
    } else if (hasClass(x,'LuHf')){
        tt <- LuHf.age(x,exterr=FALSE,i2i=i2i)
    } else if (hasClass(x,'ThU')){
        tt <- ThU.age(x,exterr=FALSE,i2i=i2i)
    }
    radialplot.default(tt,from=from,to=to,t0=t0,
                       transformation=transformation,
                       show.numbers=show.numbers,pch=pch,
                       levels=levels,clabel=clabel,bg=bg,
                       markers=markers,alpha=alpha,...)
}

radial.plot <- function(x,zeta=0,rhoD=0,asprat=3/4,
                        show.numbers=FALSE,levels=NA,clabel="",
                        bg=c('white','red'),markers=NULL,...){
    exM <- radial.scale(x,zeta,rhoD)
    tticks <- get.radial.tticks(x)
    labelpos <- 4
    if (validLevels(levels)) labelpos <- 2
    plot_radial_lines(tticks,l=0.025,x,exM[1],exM[2],
                      zeta,rhoD,label=TRUE,pos=labelpos)
    if (!is.null(markers)){
        plot_radial_lines(markers,x,exM[1],exM[2],
                          zeta,rhoD,label=FALSE)
    }
    plot_radial_axes(x)
    ns <- length(x$z)
    cols <- set.ellipse.colours(ns=ns,levels=levels,col=bg)
    plot_radial_points(x,show.numbers=show.numbers,bg=cols,...)
    colourbar(z=levels,col=bg,clabel=clabel)
}

plot_radial_points <- function(x,show.numbers=FALSE,bg='white',...){
    rxy <- data2rxry(x)
    rx <- rxy[,1]
    ry <- rxy[,2]
    if (show.numbers) {
        graphics::text(rx,ry,1:length(rx),...)
    } else if ('pch' %in% names(list(...))){
        graphics::points(rx,ry,bg=bg,...)
    } else {
        graphics::points(rx,ry,bg=bg,pch=21,...)
    }
}

plot_radial_axes <- function(x){
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
    graphics::plot(rxy[,1],rxy[,2],type='l',xlim=c(0,xM),
                   axes=FALSE,bty='n',xlab=x$xlab,
                   ylab='standardised estimate')
    c(e,xM)
}

plot_radial_lines <- function(tt,x,e,xM,zeta=0,rhoD=0,l=1,
                              label=FALSE,pos=4){
    z <- t2z(tt,x,zeta,rhoD)
    rxyb <- z2rxy(z-x$z0,e,xM)
    rxye <- z2rxy(z-x$z0,e,(1-l)*xM)
    for (i in 1:length(tt)){
        graphics::lines(c(rxyb[i,1],rxye[i,1]),
                        c(rxyb[i,2],rxye[i,2]))
        if (label) {
            if (pos==2)
                graphics::text(rxye[i,1],rxye[i,2],
                               labels=tt[i],pos=pos,xpd=NA)
            else
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

x2zs <- function(x,...){ UseMethod("x2zs",x) }
x2zs.default <- function(x,t0=NA,from=NA,to=NA,transformation=NA,...){
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
        min.z <- get.min.z(out)
        if (identical(transformation,'log')){
            out$from <- exp(min.z)-out$offset
        } else {
            out$from <- min.z
        }
    } else {
        out$from <- from
    }
    if (is.na(to)){
        max.z <- get.max.z(out)
        if (identical(transformation,'log'))
            out$to <- exp(max.z)-out$offset
        else if (identical(transformation,'linear'))
            out$to <- max.z
    } else {
        out$to <- to
    }
    out
}
x2zs.fissiontracks <- function(x,t0=NA,from=NA,to=NA,transformation=NA,...){
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
        out <- x2zs.default(tt,t0=t0,from=from,to=to,
                            transformation=transformation)
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

# this would be much easier in unicode but that doesn't render in PDF:
radial.title <- function(fit,sigdig=2){
    rounded.age <- roundit(fit$age[1],fit$age[2:3],sigdig=sigdig)
    line1 <- substitute('central age ='~a%+-%b~'|'~c,
                        list(a=rounded.age[1],
                             b=rounded.age[2],
                             c=rounded.age[3]))
    line2 <- substitute('MSWD ='~a~', p('~chi^2*')='~b,
                        list(a=signif(fit$mswd,sigdig),
                             b=signif(fit$p.value,sigdig)))
    line3 <- substitute('dispersion ='~a~'|'~b~'%',
                        list(a=signif(100*fit$disp['s'],sigdig),
                             b=signif(100*fit$disp['ci'],sigdig)))
    graphics::mtext(line1,line=2)
    graphics::mtext(line2,line=1)
    graphics::mtext(line3,line=0)
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

get.min.z <- function(x){
    rxry <- data2rxry(x)
    min.ry <- min(rxry[,2],na.rm=TRUE)
    if (min.ry > -2){ # if the data are underdispersed
        max.rx <- max(rxry[,1],na.rm=TRUE)
        min.z <- x$z0-2/max.rx
    } else {
        min.z <- min(x$z,na.rm=TRUE)
    }
    min.z
}

get.max.z <- function(x){
    rxry <- data2rxry(x)
    max.ry <- max(rxry[,2],na.rm=TRUE)
    if (max.ry < 2){ # if the data are underdispersed
        max.rx <- max(rxry[,1],na.rm=TRUE)
        max.z <- x$z0+2/max.rx
    } else {
        max.z <- max(x$z,na.rm=TRUE)
    }
    max.z
}
