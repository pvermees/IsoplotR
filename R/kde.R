#' @title
#' Create (a) kernel density estimate(s)
#' 
#' @description
#' Creates one or more kernel density estimates using a combination of
#' the Botev (2010) bandwidth selector and the Abramson (1982)
#' adaptive kernel bandwidth modifier.
#' 
#' @details
#' Given a set of \eqn{n} age estimates \eqn{\{t_1, t_2, ..., t_n\}},
#' histograms and KDEs are probability density estimators that display
#' age distributions by smoothing.  Histograms do this by grouping the
#' data into a number of regularly spaced bins.  Alternatively, kernel
#' density estimates (KDEs; Vermeesch, 2012) smooth data by applying a
#' (Gaussian) kernel:
#'
#' \eqn{KDE(t) = \sum_{i=1}^{n}N(t|\mu=t_i,\sigma=h[t])/n}
#'
#' where \eqn{N(t|\mu,\sigma)} is the probability of observing a value
#' \eqn{t} under a Normal distribution with mean \eqn{\mu} and
#' standard deviation \eqn{\sigma}.  \eqn{h[t]} is the smoothing
#' parameter or `bandwidth' of the kernel density estimate, which may
#' or may not depend on the age \eqn{t}. If \eqn{h[t]} depends on
#' \eqn{t}, then \eqn{KDE(t)} is known as an `adaptive' KDE.  The
#' default bandwidth used by \code{IsoplotR} is calculated using the
#' algorithm of Botev et al. (2010) and modulated by the adaptive
#' smoothing approach of Abramson (1982).  The rationale behind
#' adaptive kernel density estimation is to use a narrower bandwidth
#' near the peaks of the sampling distribution (where the ordered
#' dates are closely spaced in time), and a wider bandwidth in the
#' distribution's sparsely sampled troughs. Thus, the resolution of
#' the density estimate is optimised according to data availability.
#' @param x a vector of numbers OR an object of class \code{UPb},
#'     \code{PbPb}, \code{ThPb}, \code{ArAr}, \code{KCa}, \code{ReOs},
#'     \code{SmNd}, \code{RbSr}, \code{UThHe}, \code{fissiontracks},
#'     \code{ThU} or \code{detrital}
#' @rdname kde
#' @export
kde <- function(x,...){ UseMethod("kde",x) }
#' @param from minimum age of the time axis. If \code{NULL}, this is
#'     set automatically
#' @param to maximum age of the time axis. If \code{NULL}, this is set
#'     automatically
#' @param bw the bandwidth of the KDE. If \code{NULL}, \code{bw} will
#'     be calculated automatically using the algorithm by Botev et
#'     al. (2010).
#' @param adaptive logical flag controlling if the adaptive KDE
#'     modifier of Abramson (1982) is used
#' @param log transform the ages to a log scale if \code{TRUE}
#' @param n horizontal resolution (i.e., the number of segments) of
#'     the density estimate.
#' @param plot show the KDE as a plot
#' @param rug add a rug plot?
#' @param xlab the x-axis label
#' @param ylab the y-axis label
#' @param kde.col the fill colour of the KDE specified as a four
#'     element vector of \code{r, g, b, alpha} values
#' @param show.hist logical flag indicating whether a histogram should
#'     be added to the KDE
#' @param hist.col the fill colour of the histogram specified as a
#'     four element vector of \code{r, g, b, alpha} values
#' @param binwidth scalar width of the histogram bins, in Myr if
#'     \code{log = FALSE}, or as a fractional value if \code{log =
#'     TRUE}. Sturges' Rule (\eqn{\log_2[n]+1}, where \eqn{n} is the
#'     number of data points) is used if \code{binwidth = NA}
#' @param bty change to \code{"o"}, \code{"l"}, \code{"7"},
#'     \code{"c"}, \code{"u"}, or \code{"]"} if you want to draw a box
#'     around the plot
#' @param ncol scalar value indicating the number of columns over
#'     which the KDEs should be divided.
#' @param hide vector with indices of aliquots that should be removed
#'     from the plot.
#' @param nmodes label the \code{nmodes} most prominent modes of the
#'     distribution. Change to \code{'all'} to label all the modes.
#' @param sigdig the number of significant digits to which the modes
#'     should be labelled. Only used if \code{nmodes} is a positive
#'     integer or \code{'all'}.
#' @param ... optional arguments to be passed on to \code{R}'s
#'     \code{density} function.
#' @seealso \code{\link{radialplot}}, \code{\link{cad}}
#' @return If \code{x} has class \code{UPb}, \code{PbPb}, \code{ArAr},
#'     \code{KCa}, \code{ReOs}, \code{SmNd}, \code{RbSr},
#'     \code{UThHe}, \code{fissiontracks} or \code{ThU}, returns an
#'     object of class \code{KDE}, i.e. a list containing the
#'     following items:
#'
#' \describe{
#' \item{x}{ horizontal plot coordinates}
#' \item{y}{ vertical plot coordinates}
#' \item{bw}{ the base bandwidth of the density estimate}
#' \item{ages}{ the data values from the input to
#' the \code{kde} function}
#' \item{log}{ copied from the input}
#' \item{modes}{a two-column matrix with the \code{x} and \code{y}
#' values of the \code{nmodes} most prominent modes. Only returned
#' if \code{nmodes} is a positive integer or \code{'all'}.}
#' }
#'
#' or, if \code{x} has class \code{=detritals}, an object of class
#' \code{KDEs}, i.e. a list containing the following items:
#'
#' \describe{
#' \item{kdes}{a named list with objects of class \code{KDE}}
#' \item{from}{the beginning of the common time scale}
#' \item{to}{the end of the common time scale}
#' \item{themax}{the maximum probability density of all the KDEs}
#' \item{xlabel}{the x-axis label to be used by \code{plot.KDEs(...)}}
#' }
#' @references
#' Abramson, I.S., 1982. On bandwidth variation in kernel estimates-a
#' square root law. The annals of Statistics, pp.1217-1223.
#'
#' Botev, Z. I., J. F. Grotowski, and
#' D. P. Kroese. "Kernel density estimation via diffusion." The Annals
#' of Statistics 38.5 (2010): 2916-2957.
#'
#' Vermeesch, P., 2012. On the visualisation of detrital age
#' distributions. Chemical Geology, 312, pp.190-194.
#' @examples
#' kde(examples$UPb)
#'
#' dev.new()
#' kde(examples$FT1,log=TRUE)
#'
#' dev.new()
#' kde(examples$DZ,from=1,to=3000,kernel="epanechnikov")
#' @importFrom grDevices rgb
#' @rdname kde
#' @export
kde.default <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,log=FALSE,
                        n=512,plot=TRUE,rug=TRUE,xlab="age [Ma]",
                        ylab="",kde.col=rgb(1,0,1,0.6),
                        hist.col=rgb(0,1,0,0.2),show.hist=TRUE,bty='n',
                        binwidth=NA,hide=NULL,nmodes=0,sigdig=2,...){
    x2calc <- clear(x,hide)
    X <- get_kde(x2calc,from=from,to=to,bw=bw,adaptive=adaptive,
                log=log,n=n,nmodes=nmodes,...)
    if (plot) {
        plot.KDE(X,rug=rug,xlab=xlab,ylab=ylab,kde.col=kde.col,
                 hist.col=hist.col,show.hist=show.hist,bty=bty,
                 binwidth=binwidth,sigdig=sigdig)
    }
    invisible(X)
}
#' @rdname kde
#' @export
kde.other <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,log=FALSE,
                      n=512,plot=TRUE,rug=TRUE,xlab="age [Ma]",
                      ylab="",kde.col=rgb(1,0,1,0.6),
                      hist.col=rgb(0,1,0,0.2),show.hist=TRUE,
                      bty='n',binwidth=NA,hide=NULL,
                      nmodes=0,sigdig=2,...){
    if (x$format<3) X <- x$x[,1]
    else if (x$format==3) X <- x$x[,2]
    else stop("KDEs are not available for this format")
    kde(X,from=from,to=to,bw=bw,adaptive=adaptive,log=log,
        n=n,plot=plot,rug=rug,xlab=xlab,ylab=ylab,kde.col=kde.col,
        hist.col=hist.col,show.hist=show.hist,bty=bty,
        binwidth=binwidth,hide=hide,nmodes=nmodes,sigdig=sigdig,...)
}
#' @param type scalar indicating whether to plot the
#'     \eqn{^{207}}Pb/\eqn{^{235}}U age (\code{type}=1), the
#'     \eqn{^{206}}Pb/\eqn{^{238}}U age (\code{type}=2), the
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb age (\code{type}=3), the
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb-\eqn{^{206}}Pb/\eqn{^{238}}U age
#'     (\code{type}=4), the concordia_age (\code{type}=5), or the
#'     \eqn{^{208}}Pb/\eqn{^{232}}Th age (\code{type}=6).
#' @param cutoff.76 the age (in Ma) below which the
#'     \eqn{^{206}}Pb/\eqn{^{238}}U and above which the
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb age is used. This parameter is
#'     only used if \code{type=4}.
#' @param cutoff.disc discordance cutoff filter. This is an object of
#'     class \code{\link{discfilter}}.
#' 
#' @param common.Pb common lead correction:
#'
#' \code{0}: none
#'
#' \code{1}: use the Pb-composition stored in
#' 
#' \code{settings('iratio','Pb207Pb206')} (if \code{x} has class
#' \code{UPb} and \code{x$format<4});
#' 
#' \code{settings('iratio','Pb206Pb204')} and
#' \code{settings('iratio','Pb207Pb204')} (if \code{x} has class
#' \code{PbPb} or \code{x} has class \code{UPb} and
#' \code{3<x$format<7}); or
#'
#' \code{settings('iratio','Pb206Pb208')} and
#' \code{settings('iratio','Pb207Pb208')} (if \code{x} has class
#' \code{UPb} and \code{x$format=7,8}).
#'
#' \code{2}: use the isochron intercept as the initial Pb-composition
#'
#' \code{3}: use the Stacey-Kramers two-stage model to infer the
#' initial Pb-composition (only valid if \code{x} has class
#' \code{UPb}).
#' @rdname kde
#' @export
kde.UPb <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,log=FALSE,
                    n=512,plot=TRUE,rug=TRUE,xlab="age [Ma]",ylab="",
                    kde.col=rgb(1,0,1,0.6),hist.col=rgb(0,1,0,0.2),
                    show.hist=TRUE, bty='n',binwidth=NA,type=4,
                    cutoff.76=1100,cutoff.disc=discfilter(),
                    common.Pb=0,hide=NULL,nmodes=0,sigdig=2,...){
    kde_helper(x,type=type,cutoff.76=cutoff.76,
               cutoff.disc=cutoff.disc,common.Pb=common.Pb,
               from=from,to=to,bw=bw,adaptive=adaptive,log=log,
               n=n,plot=plot,rug=rug,xlab=xlab,ylab=ylab,
               kde.col=kde.col,hist.col=hist.col,
               show.hist=show.hist,bty=bty,binwidth=binwidth,
               hide=hide,nmodes=nmodes,sigdig=sigdig,...)
}
#' @param samebandwidth logical flag indicating whether the same
#'     bandwidth should be used for all samples. If
#'     \code{samebandwidth = TRUE} and \code{bw = NULL}, then the
#'     function will use the median bandwidth of all the samples.
#' 
#' @param normalise logical flag indicating whether or not the KDEs
#'     should all integrate to the same value.
#' 
#' @rdname kde
#' @export
kde.detritals <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,
                          log=FALSE, n=512,plot=TRUE,rug=FALSE,
                          xlab="age [Ma]",ylab="",
                          kde.col=rgb(1,0,1,0.6),
                          hist.col=rgb(0,1,0,0.2),show.hist=TRUE,
                          bty='n',binwidth=NA,ncol=NA,
                          samebandwidth=TRUE,normalise=TRUE,
                          hide=NULL,nmodes=0,sigdig=2,...){
    if (is.character(hide)) hide <- which(names(x)%in%hide)
    x2plot <- clear(x,hide)
    X <- get_kde(x2plot,from=from,to=to,bw=bw,adaptive=adaptive,
                log=log,n=n,samebandwidth=samebandwidth,
                normalise=normalise,nmodes=nmodes,...)
    if (plot){
        plot.KDEs(X,rug=rug,xlab=xlab,ylab=ylab,kde.col=kde.col,
                  hist.col=hist.col,show.hist=show.hist,bty=bty,
                  binwidth=binwidth,ncol=ncol,sigdig=sigdig)
    }
    invisible(X)
}
#' @param i2i `isochron to intercept': calculates the initial (aka
#'     `inherited', `excess', or `common')
#'     \eqn{^{40}}Ar/\eqn{^{36}}Ar, \eqn{^{40}}Ca/\eqn{^{44}}Ca,
#'     \eqn{^{207}}Pb/\eqn{^{204}}Pb, \eqn{^{87}}Sr/\eqn{^{86}}Sr,
#'     \eqn{^{143}}Nd/\eqn{^{144}}Nd, \eqn{^{187}}Os/\eqn{^{188}}Os,
#'     \eqn{^{230}}Th/\eqn{^{232}}Th, \eqn{^{176}}Hf/\eqn{^{177}}Hf or
#'     \eqn{^{204}}Pb/\eqn{^{208}}Pb ratio from an isochron
#'     fit. Setting \code{i2i} to \code{FALSE} uses the default values
#'     stored in \code{settings('iratio',...)}.
#' @rdname kde
#' @export
kde.PbPb <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,log=FALSE,
                     n=512,plot=TRUE,rug=TRUE,xlab="age [Ma]",ylab="",
                     kde.col=rgb(1,0,1,0.6),hist.col=rgb(0,1,0,0.2),
                     show.hist=TRUE,bty='n',binwidth=NA,common.Pb=2,
                     hide=NULL,nmodes=0,sigdig=2,...){
    kde_helper(x,common.Pb=common.Pb,from=from,to=to,bw=bw,
               adaptive=adaptive,log=log,n=n,plot=plot,rug=rug,
               xlab=xlab,ylab=ylab,kde.col=kde.col,hist.col=hist.col,
               show.hist=show.hist,bty=bty,binwidth=binwidth,
               hide=hide,nmodes=nmodes,sigdig=sigdig,...)
}
#' @rdname kde
#' @export
kde.ArAr <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,log=FALSE,
                     n=512,plot=TRUE,rug=TRUE,xlab="age [Ma]",ylab="",
                     kde.col=rgb(1,0,1,0.6),hist.col=rgb(0,1,0,0.2),
                     show.hist=TRUE,bty='n',binwidth=NA,i2i=FALSE,
                     hide=NULL,nmodes=0,sigdig=2,...){
    kde_helper(x,i2i=i2i,from=from,to=to,bw=bw,adaptive=adaptive,
               log=log,n=n,plot=plot,rug=rug,xlab=xlab,ylab=ylab,
               kde.col=kde.col,hist.col=hist.col,show.hist=show.hist,
               bty=bty,binwidth=binwidth,hide=hide,
               nmodes=nmodes,sigdig=sigdig,...)
}
#' @rdname kde
#' @export
kde.KCa <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,log=FALSE,
                    n=512,plot=TRUE,rug=TRUE,xlab="age [Ma]",ylab="",
                    kde.col=rgb(1,0,1,0.6),hist.col=rgb(0,1,0,0.2),
                    show.hist=TRUE,bty='n',binwidth=NA,i2i=FALSE,
                    hide=NULL,nmodes=0,sigdig=2,...){
    kde_helper(x,i2i=i2i,from=from,to=to,bw=bw,adaptive=adaptive,
               log=log,n=n,plot=plot,rug=rug,xlab=xlab,ylab=ylab,
               kde.col=kde.col,hist.col=hist.col,show.hist=show.hist,
               bty=bty,binwidth=binwidth,hide=hide,
               nmodes=nmodes,sigdig=sigdig,...)
}
#' @rdname kde
#' @export
kde.ThPb <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,log=FALSE,
                     n=512,plot=TRUE,rug=TRUE,xlab="age [Ma]",ylab="",
                     kde.col=rgb(1,0,1,0.6),hist.col=rgb(0,1,0,0.2),
                     show.hist=TRUE,bty='n',binwidth=NA,i2i=FALSE,
                     hide=NULL,nmodes=0,sigdig=2,...){
    kde_helper(x,i2i=i2i,from=from,to=to,bw=bw,adaptive=adaptive,
               log=log,n=n,plot=plot,rug=rug,xlab=xlab,ylab=ylab,
               kde.col=kde.col,hist.col=hist.col,show.hist=show.hist,
               bty=bty,binwidth=binwidth,hide=hide,
               nmodes=nmodes,sigdig=sigdig,...)
}
#' @param Th0i initial \eqn{^{230}}Th correction.
#'
#' \code{0}: no correction
#'
#' \code{1}: project the data along an isochron fit
#'
#' \code{2}: if \code{x$format} is \code{1} or \code{2}, correct the
#' data using the measured present day \eqn{^{230}}Th/\eqn{^{238}}U,
#' \eqn{^{232}}Th/\eqn{^{238}}U and \eqn{^{234}}U/\eqn{^{238}}U
#' activity ratios in the detritus. If \code{x$format} is \code{3} or
#' \code{4}, correct the data using the measured
#' \eqn{^{238}}U/\eqn{^{232}}Th activity ratio of the whole rock, as
#' stored in \code{x} by the \code{read.data()} function.
#'
#' \code{3}: correct the data using an assumed initial
#' \eqn{^{230}}Th/\eqn{^{232}}Th-ratio for the detritus (only relevant
#' if \code{x$format} is \code{1} or \code{2}).
#'
#' @rdname kde
#' @export
kde.ThU <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,log=FALSE,
                    n=512,plot=TRUE,rug=TRUE,xlab="age [ka]",ylab="",
                    kde.col=rgb(1,0,1,0.6),hist.col=rgb(0,1,0,0.2),
                    show.hist=TRUE,bty='n',binwidth=NA,Th0i=0,hide=NULL,
                    nmodes=0,sigdig=2,...){
    kde_helper(x,Th0i=Th0i,from=from,to=to,bw=bw,adaptive=adaptive,
               log=log,n=n,plot=plot,rug=rug,xlab=xlab,ylab=ylab,
               kde.col=kde.col,hist.col=hist.col,show.hist=show.hist,
               bty=bty,binwidth=binwidth,hide=hide,
               nmodes=nmodes,sigdig=sigdig,...)
}
#' @rdname kde
#' @export
kde.ReOs <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,log=FALSE,
                     n=512,plot=TRUE,rug=TRUE,xlab="age [Ma]",ylab="",
                     kde.col=rgb(1,0,1,0.6),hist.col=rgb(0,1,0,0.2),
                     show.hist=TRUE,bty='n',binwidth=NA,i2i=TRUE,
                     hide=NULL,nmodes=0,sigdig=2,...){
    kde_helper(x,i2i=i2i,from=from,to=to,bw=bw,adaptive=adaptive,
               log=log,n=n,plot=plot,rug=rug,xlab=xlab,ylab=ylab,
               kde.col=kde.col,hist.col=hist.col,show.hist=show.hist,
               bty=bty,binwidth=binwidth,hide=hide,
               nmodes=nmodes,sigdig=sigdig,...)
}
#' @rdname kde
#' @export
kde.SmNd <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,log=FALSE,
                     n=512,plot=TRUE,rug=TRUE,xlab="age [Ma]",ylab="",
                     kde.col=rgb(1,0,1,0.6),hist.col=rgb(0,1,0,0.2),
                     show.hist=TRUE,bty='n',binwidth=NA,i2i=TRUE,
                     hide=NULL,nmodes=0,sigdig=2,...){
    kde_helper(x,i2i=i2i,from=from,to=to,bw=bw,adaptive=adaptive,
               log=log,n=n,plot=plot,rug=rug,xlab=xlab,ylab=ylab,
               kde.col=kde.col,hist.col=hist.col,show.hist=show.hist,
               bty=bty,binwidth=binwidth,hide=hide,
               nmodes=nmodes,sigdig=sigdig,...)
}
#' @rdname kde
#' @export
kde.RbSr <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,log=FALSE,
                     n=512,plot=TRUE,rug=TRUE,xlab="age [Ma]",ylab="",
                     kde.col=rgb(1,0,1,0.6),hist.col=rgb(0,1,0,0.2),
                     show.hist=TRUE,bty='n',binwidth=NA,i2i=TRUE,
                     hide=NULL,nmodes=0,sigdig=2,...){
    kde_helper(x,i2i=i2i,from=from,to=to,bw=bw,adaptive=adaptive,
               log=log,n=n,plot=plot,rug=rug,xlab=xlab,ylab=ylab,
               kde.col=kde.col,hist.col=hist.col,show.hist=show.hist,
               bty=bty,binwidth=binwidth,hide=hide,
               nmodes=nmodes,sigdig=sigdig,...)
}
#' @rdname kde
#' @export
kde.LuHf <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,log=FALSE,
                     n=512,plot=TRUE,rug=TRUE,xlab="age [Ma]",ylab="",
                     kde.col=rgb(1,0,1,0.6),hist.col=rgb(0,1,0,0.2),
                     show.hist=TRUE,bty='n',binwidth=NA,i2i=TRUE,
                     hide=NULL,nmodes=0,sigdig=2,...){
    kde_helper(x,i2i=i2i,from=from,to=to,bw=bw,adaptive=adaptive,
               log=log,n=n,plot=plot,rug=rug,xlab=xlab,ylab=ylab,
               kde.col=kde.col,hist.col=hist.col,show.hist=show.hist,
               bty=bty,binwidth=binwidth,hide=hide,
               nmodes=nmodes,sigdig=sigdig,...)
}
#' @rdname kde
#' @export
kde.UThHe <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,log=FALSE,
                      n=512,plot=TRUE,rug=TRUE,xlab="age [Ma]",ylab="",
                      kde.col=rgb(1,0,1,0.6),hist.col=rgb(0,1,0,0.2),
                      show.hist=TRUE,bty='n',binwidth=NA,
                      hide=NULL,nmodes=0,sigdig=2,...){
    kde_helper(x,from=from,to=to,bw=bw,adaptive=adaptive,log=log,
               n=n,plot=plot,rug=rug,xlab=xlab,ylab=ylab,
               kde.col=kde.col, hist.col=hist.col,show.hist=show.hist,
               bty=bty,binwidth=binwidth,hide=hide,
               nmodes=nmodes,sigdig=sigdig,...)
}
#' @rdname kde
#' @export
kde.fissiontracks <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,
                              log=FALSE,n=512,plot=TRUE,rug=TRUE,
                              xlab="age [Ma]",ylab="",
                              kde.col=rgb(1,0,1,0.6),
                              hist.col=rgb(0,1,0,0.2),show.hist=TRUE,
                              bty='n',binwidth=NA,hide=NULL,
                              nmodes=0,sigdig=2,...){
    kde_helper(x,from=from,to=to,bw=bw,adaptive=adaptive,log=log,
               n=n,plot=plot,rug=rug,xlab=xlab,ylab=ylab,
               kde.col=kde.col,hist.col=hist.col,show.hist=show.hist,
               bty=bty,binwidth=binwidth,hide=hide,
               nmodes=nmodes,sigdig=sigdig,...)
}

kde_helper <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,log=FALSE,
                       n=512,plot=TRUE,rug=TRUE,xlab="age [Ma]",ylab="",
                       kde.col=rgb(1,0,1,0.6),hist.col=rgb(0,1,0,0.2),
                       show.hist=TRUE, bty='n',binwidth=NA,type=4,
                       cutoff.76=1100,cutoff.disc=discfilter(),
                       common.Pb=0,i2i=FALSE,Th0i=0,hide=NULL,
                       nmodes=0,sigdig=2,...){
    tt <- get_ages(x,type=type,cutoff.76=cutoff.76,
                   cutoff.disc=cutoff.disc,i2i=i2i,
                   common.Pb=common.Pb,Th0i=Th0i,omit4c=hide)
    kde.default(tt[,1],from=from,to=to,bw=bw,adaptive=adaptive,
                log=log,n=n,plot=plot,rug=rug,xlab=xlab,ylab=ylab,
                kde.col=kde.col,hist.col=hist.col,hide=hide,
                show.hist=show.hist,bty=bty,binwidth=binwidth,
                nmodes=nmodes,sigdig=sigdig,...)
}

#' Helper functions for the generic kde function
#' @param x an IsoplotR object or a vector of numbers
#' @noRd
get_kde <- function(x,...){ UseMethod("get_kde",x) }
#' @param from minimum limit
#' @param to maximum limit
#' @param bw bandwidth
#' @param adaptive if \code{TRUE}, adjusts the bandwidth by density
#' @param log if \code{TRUE}, forms a KDE of the logarithms
#' @param n the number of steps at which the KDE is evaluated
#' @param nmodes the number of modes to be labelled
#' @noRd
get_kde.default <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,
                           log=FALSE,n=512,nmodes=0,...){
    out <- list()
    class(out) <- "KDE"
    out$name <- deparse(substitute(x))
    out$log <- log
    if (log) d <- log(x)
    else d <- x
    if (is.na(bw)) bw <- botev(d)
    if (is.na(from) | is.na(to)) {
        mM <- getmM(x,from,to,log)
        to <- mM$M + bw
        from <- mM$m
        if (mM$m > bw) from <- from - bw
    }
    if (log) {
        from <- log(from)
        to <- log(to)
    }
    xx <- seq(from=from,to=to,length.out=n)
    if (adaptive){
        yy <- Abramson(d,from=from,to=to,bw=bw,n=n,...)
    } else {
        yy <- stats::density(d,bw,from=from,to=to,n=n,...)$y
    }
    if (log){ # back-transform
        dxin <- diff(xx)
        out$x <- exp(xx)
        dxout <- diff(out$x)
        yy <- c(dxin,utils::tail(dxin,1))*yy/c(dxout,utils::tail(dxout,1))
    } else {
        out$x <- xx
    }
    out$x <- c(out$x[1],out$x,out$x[n])
    out$y <- c(0,yy/(sum(yy)*(to-from)/n),0)
    out$modes <- getmodes(x=out$x,y=out$y,nmodes=nmodes)
    out$bw <- bw
    out$ages <- x
    out
}
#' @param samebandwidth use the median bandwidth for all samples?
#' @param normalise rescale the KDEs to cover the same area
#' @noRd
get_kde.detritals <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,
                             log=FALSE,n=512,samebandwidth=TRUE,
                             normalise=TRUE,nmodes=0,...){
    if (is.na(from) | is.na(to)) {
        mM <- getmM(unlist(x),from,to,log)
        from <- mM$m
        to <- mM$M
    }
    snames <- names(x)
    thekdes <- list()
    themax <- -1
    if (is.na(bw) & samebandwidth) bw <- commonbandwidth(x,log=log)
    for (name in snames){
        thekdes[[name]] <- kde(x[[name]],from=from,to=to,bw=bw,
                               adaptive=adaptive,log=log,n=n,
                               plot=FALSE,nmodes=nmodes,...)
        if (normalise){
            maxval <- max(thekdes[[name]]$y)
            if (themax < maxval) {themax <- maxval}
        }
    }
    out <- list()
    out$kdes <- thekdes
    out$from <- from
    out$to <- to
    out$themax <- themax
    out$log <- log
    class(out) <- "KDEs"
    out
}

getmodes <- function(x,y,miny=max(y)/1000,nmodes=0){
    if (nmodes==0){
        return(NULL)
    } else {
        slope <- diff(y)
        shape <- diff(sign(slope))
        imax <- which(shape<0)+1
        imodes <- imax[order(y[imax],decreasing=TRUE)]
        toosmall <- y[imodes] < miny
        modes <- cbind(x=x[imodes[!toosmall]],y=y[imodes[!toosmall]])
    }
    if (identical(nmodes,'all')){
        out <- modes
    } else {
        nm <- min(c(nrow(modes),nmodes))
        out <- modes[1:nm,,drop=FALSE]
    }
    out
}

# Abramson
# get geometric mean pilot density
getG <- function(pdens) {
    fpos <- pdens[pdens>0]
    N <- length(fpos)
    logmean <- mean(log(fpos))
    exp(logmean)
}

# Abramson
# get fixed bandwidth pilot density
pilotdensity <- function(dat,bw){
    d <- stats::na.omit(dat)
    n <- length(d)
    dens <- rep(0,n)
    for (i in 1:n){
        dens <- dens + stats::dnorm(d,mean=d[i],sd=bw)
    }
    dens
}

# adaptive KDE algorithm of Abramson (1982) as summarised by Jahn (2007)
Abramson <- function(dat,from,to,bw,n=512,...){
    d <- stats::na.omit(dat)
    nn <- length(d)
    pdens <- pilotdensity(d,bw)
    G <- getG(pdens)
    lambda <- 0
    dens <- rep(0,n)
    for (i in 1:nn){
        lambda = sqrt(G/pdens[i])
        dens <- dens + stats::density(d[i],bw*lambda,
                                      from=from,to=to,n=n,...)$y
    }
    dens
}

#' @noRd
plot.KDE <- function(x,rug=TRUE,xlab="age [Ma]",ylab="",
                     kde.col=rgb(1,0,1,0.6),show.hist=TRUE,
                     hist.col=rgb(0,1,0,0.2),
                     binwidth=NA,bty='n',sigdig=2,...){
    m <- x$x[1]
    M <- utils::tail(x$x,n=1)
    inrange <- x$ages >= m & x$ages <= M
    ages <- x$ages[inrange]
    if (is.na(binwidth)) nb <- log2(length(ages))+1 # Sturges' Rule
    if (x$log){
        do.log <- 'x'
        if (is.na(binwidth)) {
            breaks <- exp(seq(log(m),log(M),length.out=nb+1))
        } else {
            breaks <- exp(seq(log(m),log(M)+binwidth,by=binwidth))
        }
        if (M/m < breaks[2]/breaks[1]) show.hist <- FALSE
        else h <- graphics::hist(log(ages),breaks=log(breaks),plot=FALSE)
    } else {
        do.log <- ''
        if (is.na(binwidth)) {
            breaks <- seq(m,M,length.out=nb+1)
        } else {
            breaks <- seq(m,M+binwidth,by=binwidth)
        }
        h <- graphics::hist(ages,breaks=breaks,plot=FALSE)
    }
    nb <- length(breaks)-1
    maxy <- max(x$y)
    if (show.hist) maxy <- max(maxy,max(h$density))
    graphics::plot(range(x$x),c(0,maxy),type='n',log=do.log,
                   xlab=xlab,ylab=ylab,yaxt='n',bty=bty,...)
    if (show.hist){
        graphics::rect(xleft=breaks[1:nb],xright=breaks[2:(nb+1)],
                       ybottom=0,ytop=h$density,col=hist.col)
        if (graphics::par('yaxt')!='n') {
            fact <- max(h$counts)/max(h$density)
            lbls <- pretty(fact*h$density)
            at <- lbls/fact
            graphics::axis(2,at=at,labels=lbls)
        }
    }
    graphics::polygon(x$x,x$y,col=kde.col)
    graphics::lines(x$x,x$y,col='black')
    if (rug) rug(ages)
    if (!is.null(x$modes)){
        graphics::points(x=x$modes[,'x'],y=x$modes[,'y'],pch='|')
        graphics::text(x=x$modes[,'x'],y=x$modes[,'y'],xpd=TRUE,pos=3,
                       labels=roundit(x$modes[,'x'],sigdig=sigdig))
    }
    mymtext(get_ntit(x$ages,m=m,M=M),line=0,adj=1)
}

#' @noRd
plot.KDEs <- function(x,ncol=NA,rug=FALSE,xlab="age [Ma]",ylab="",
                      kde.col=rgb(1,0,1,0.6),show.hist=TRUE,
                      hist.col=rgb(0,1,0,0.2),
                      binwidth=NA,bty='n',sigdig=2,...){
    if (is.na(ncol)) ncol <- ceiling(sqrt(length(x)/2))
    oldpar <- graphics::par(no.readonly=T)
    snames <- names(x$kdes)
    ns <- length(snames)
    w <- rep(1,ncol) # column widths
    nppc <- ceiling(ns/ncol)
    np <- nppc*ncol # number of subpanels
    graphics::layout(matrix(1:np,nppc,length(w)),w,rep(1,nppc))
    si <- ceiling(seq(from=0,to=ns,length.out=ncol+1)) # sample index
    graphics::par(xpd=TRUE,mar=rep(1,4),oma=c(3,1,1,1))
    if (x$themax>0) ylim <- c(0,x$themax)
    else ylim <- NULL
    for (i in 1:ns){
        if ((i%%nppc)==0 | (i==ns)) {
            plot.KDE(x$kdes[[i]],rug=rug,xlab=xlab,ylab=ylab,
                     kde.col=kde.col,show.hist=show.hist,
                     hist.col=hist.col,binwidth=binwidth,
                     bty=bty,ann=FALSE,ylim=ylim,sigdig=sigdig,...)
            graphics::mtext(side=1,text=xlab,line=2)
        } else {
            plot.KDE(x$kdes[[i]],rug=rug,xlab=xlab,ylab=ylab,
                     kde.col=kde.col,show.hist=show.hist,
                     hist.col=hist.col,binwidth=binwidth,
                     bty=bty,xaxt='n',ylim=ylim,sigdig=sigdig,...)
        }
        graphics::title(snames[i])
    }
    graphics::par(oldpar)
}
