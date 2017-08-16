#' Create (a) kernel density estimate(s)
#'
#' Creates one or more kernel density estimates using a combination of
#' the Botev (2010) bandwidth selector and the Abramson (1982)
#' adaptive kernel bandwidth modifier.
#'
#' @param x a vector of numbers OR an object of class \code{UPb},
#'     \code{PbPb}, \code{ArAr}, \code{ReOs}, \code{SmNd},
#'     \code{RbSr}, \code{UThHe}, \code{fissiontracks}, \code{ThU} or
#'     \code{detrital}
#' @rdname kde
#' @export
kde <- function(x,...){ UseMethod("kde",x) }
#' @param from minimum age of the time axis. If \code{NULL}, this is set
#'     automatically
#' @param to maximum age of the time axis. If \code{NULL}, this is set
#'     automatically
#' @param bw the bandwidth of the KDE. If \code{NULL}, \code{bw} will be calculated
#'     automatically using \code{botev()}
#' @param adaptive logical flag controlling if the adaptive KDE
#'     modifier of Abramson (1982) is used
#' @param log transform the ages to a log scale if \code{TRUE}
#' @param n horizontal resolution of the density estimate
#' @param plot show the KDE as a plot
#' @param pch the symbol used to show the samples. May be a vector.
#'     Set \code{pch=NA} to turn them off.
#' @param xlab the x-axis label
#' @param ylab the y-axis label
#' @param kde.col the fill colour of the KDE specified as a four
#'     element vector of \code{r, g, b, alpha} values
#' @param show.hist logical flag indicating whether a histogram should
#'     be added to the KDE
#' @param hist.col the fill colour of the histogram specified as a
#'     four element vector of r, g, b, alpha values
#' @param binwidth scalar width of the histogram bins, in Myr if
#'     \code{x$log = FALSE}, or as a fractional value if
#'     \code{x$log = TRUE}. Sturges' Rule is used if
#'     \code{binwidth = NA}
#' @param bty change to \code{"o"}, \code{"l"}, \code{"7"},
#'     \code{"c"}, \code{"u"}, or \code{"]"} if you want to draw a box
#'     around the plot
#' @param ncol scalar value indicating the number of columns over
#'     which the KDEs should be divided. This option is only used if
#'     \code{x} has class \code{detritals}.
#' @param ... optional arguments to be passed on to \code{density}
#' @return
#'
#' if \code{plot = TRUE}, returns an object of class \code{KDE}, i.e. a
#'     list containing the following items:
#' 
#' \describe{
#' \item{x}{ horizontal plot coordinates}
#' \item{y}{ vertical plot coordinates}
#' \item{bw}{ the base bandwidth of the density estimate}
#' \item{ages}{ the data values from the input to the \code{kde} function}
#' }
#'
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
#'
#' @examples
#' data(examples)
#' kde(examples$DZ[['N1']],kernel="epanechnikov")
#' @importFrom grDevices rgb
#' @rdname kde
#' @export
kde.default <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,log=FALSE,
                        n=512,plot=TRUE,pch=NA,xlab="age [Ma]",
                        ylab="",kde.col=rgb(1,0,1,0.6),
                        hist.col=rgb(0,1,0,0.2),show.hist=TRUE,
                        bty='n',binwidth=NA,ncol=NA,...){
    X <- getkde(x,from=from,to=to,bw=bw,adaptive=adaptive,log=log,n=n,...)
    if (plot) {
        plot.KDE(X,pch=pch,xlab=xlab,ylab=ylab,kde.col=kde.col,
                 hist.col=hist.col,show.hist=show.hist,bty=bty,
                 binwidth=binwidth)
    } else {
        return(X)
    }
}
#' @param type scalar indicating whether to plot the
#'     \eqn{^{207}}Pb/\eqn{^{235}}U age (\code{type}=1), the
#'     \eqn{^{206}}Pb/\eqn{^{238}}U age (\code{type}=2), the
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb age (\code{type}=3), the
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
#' @rdname kde
#' @export
kde.UPb <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,log=FALSE,
                    n=512,plot=TRUE,pch=NA,xlab="age [Ma]",ylab="",
                    kde.col=rgb(1,0,1,0.6),hist.col=rgb(0,1,0,0.2),
                    show.hist=TRUE, bty='n',binwidth=NA,ncol=NA,
                    type=4,cutoff.76=1100,cutoff.disc=c(-15,5),
                    common.Pb=0,...){
    if (common.Pb %in% c(1,2,3))
        X <- common.Pb.correction(x,option=common.Pb)
    else
        X <- x
    tt <- filter.UPb.ages(X,type,cutoff.76,cutoff.disc)[,1]
    kde.default(tt,from=from,to=to,bw=bw,adaptive=adaptive,log=log,
                n=n,plot=plot,pch=pch,xlab=xlab,ylab=ylab,
                kde.col=kde.col, hist.col=hist.col,
                show.hist=show.hist,bty=bty,binwidth=binwidth,
                ncol=ncol,...)
}
#' @param samebandwidth logical flag indicating whether the same
#'     bandwidth should be used for all samples. If
#'     \code{samebandwidth = TRUE} and \code{bw = NULL}, then the
#'     function will use the median bandwidth of all the samples.
#' @param normalise logical flag indicating whether or not the KDEs
#'     should all integrate to the same value.
#' @return
#'
#' or, if \code{x} has class \code{=detritals}, an object of class
#' \code{KDEs}, i.e. a list containing the following items:
#'
#' \describe{
#' \item{kdes}{a named list with objects of class \code{KDE}}
#' \item{from}{the beginning of the common time scale}
#' \item{to}{the end of the common time scale}
#' \item{themax}{the maximum probability density of all the KDEs}
#' \item{xlabel}{the x-axis label to be used by \code{plot.KDEs}}
#' }
#' 
#' @examples
#' kde(examples$DZ,from=0,to=3000)
#' @rdname kde
#' @export
kde.detritals <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,
                          log=FALSE, n=512,plot=TRUE,pch=NA,
                          xlab="age [Ma]",ylab="",
                          kde.col=rgb(1,0,1,0.6),
                          hist.col=rgb(0,1,0,0.2),show.hist=TRUE,
                          bty='n',binwidth=NA,ncol=NA,
                          samebandwidth=TRUE,normalise=TRUE,...){
    X <- getkde(x,from=from,to=to,bw=bw,adaptive=adaptive,log=log,n=n,
                samebandwidth=samebandwidth,normalise=normalise,...)
    if (plot){
        plot.KDEs(X,pch=pch,xlab=xlab,ylab=ylab,kde.col=kde.col,
                  hist.col=hist.col,show.hist=show.hist,bty=bty,
                  binwidth=binwidth,ncol=ncol)
    } else {
        return(X)
    }
}
#' @param i2i `isochron to intercept': calculates the initial
#'     (aka `inherited', `excess', or `common') \eqn{^{40}}Ar/\eqn{^{36}}Ar,
#'     \eqn{^{207}}Pb/\eqn{^{204}}Pb, \eqn{^{87}}Sr/\eqn{^{86}}Sr,
#'     \eqn{^{143}}Nd/\eqn{^{144}}Nd, \eqn{^{187}}Os/\eqn{^{188}}Os or
#'     \eqn{^{176}}Hf/\eqn{^{177}}Hf ratio from an isochron
#'     fit. Setting \code{i2i} to \code{FALSE} uses the default values
#'     stored in \code{settings('iratio',...)}  or zero (for the Pb-Pb
#'     method). When applied to data of class \code{ThU}, setting
#'     \code{i2i} to \code{TRUE} applies a detrital Th-correction.
#' @rdname kde
#' @export
kde.PbPb <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,log=FALSE,
                     n=512,plot=TRUE,pch=NA,xlab="age [Ma]",ylab="",
                     kde.col=rgb(1,0,1,0.6),hist.col=rgb(0,1,0,0.2),
                     show.hist=TRUE,bty='n',binwidth=NA,ncol=NA,
                     i2i=FALSE,...){
    tt <- PbPb.age(x,i2i=i2i)[,1]
    kde.default(tt,from=from,to=to,bw=bw,adaptive=adaptive,log=log,
                n=n,plot=plot,pch=pch,xlab=xlab,ylab=ylab,
                kde.col=kde.col,hist.col=hist.col,
                show.hist=show.hist,bty=bty,binwidth=binwidth,
                ncol=ncol,...)
}
#' @rdname kde
#' @export
kde.ArAr <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,log=FALSE,
                     n=512,plot=TRUE,pch=NA,xlab="age [Ma]",ylab="",
                     kde.col=rgb(1,0,1,0.6),hist.col=rgb(0,1,0,0.2),
                     show.hist=TRUE,bty='n',binwidth=NA,ncol=NA,
                     i2i=FALSE,...){
    tt <- ArAr.age(x,i2i=i2i)[,1]
    kde.default(tt,from=from,to=to,bw=bw,adaptive=adaptive,log=log,
                n=n,plot=plot,pch=pch,xlab=xlab,ylab=ylab,
                kde.col=kde.col,hist.col=hist.col,
                show.hist=show.hist,bty=bty,binwidth=binwidth,
                ncol=ncol,...)
}
#' @rdname kde
#' @export
kde.ThU <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,log=FALSE,
                    n=512,plot=TRUE,pch=NA,xlab="age [ka]",ylab="",
                    kde.col=rgb(1,0,1,0.6),hist.col=rgb(0,1,0,0.2),
                    show.hist=TRUE,bty='n',binwidth=NA,ncol=NA,
                    i2i=FALSE,...){
    tt <- ThU.age(x,i2i=i2i)[,1]
    kde.default(tt,from=from,to=to,bw=bw,adaptive=adaptive,log=log,
                n=n,plot=plot,pch=pch,xlab=xlab,ylab=ylab,
                kde.col=kde.col,hist.col=hist.col,
                show.hist=show.hist,bty=bty,binwidth=binwidth,
                ncol=ncol,...)
}
#' @rdname kde
#' @export
kde.ReOs <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,log=FALSE,
                     n=512,plot=TRUE,pch=NA,xlab="age [Ma]",ylab="",
                     kde.col=rgb(1,0,1,0.6),hist.col=rgb(0,1,0,0.2),
                     show.hist=TRUE,bty='n',binwidth=NA,ncol=NA,
                     i2i=TRUE,...){
    tt <- ReOs.age(x,i2i=i2i)[,1]
    kde.default(tt,from=from,to=to,bw=bw,adaptive=adaptive,log=log,
                n=n,plot=plot,pch=pch,xlab=xlab,ylab=ylab,
                kde.col=kde.col,hist.col=hist.col,
                show.hist=show.hist,bty=bty,binwidth=binwidth,
                ncol=ncol,...)
}
#' @rdname kde
#' @export
kde.SmNd <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,log=FALSE,
                     n=512,plot=TRUE,pch=NA,xlab="age [Ma]",ylab="",
                     kde.col=rgb(1,0,1,0.6),hist.col=rgb(0,1,0,0.2),
                     show.hist=TRUE,bty='n',binwidth=NA,ncol=NA,
                     i2i=TRUE,...){
    tt <- SmNd.age(x,i2i=i2i)[,1]
    kde.default(tt,from=from,to=to,bw=bw,adaptive=adaptive,log=log,
                n=n,plot=plot,pch=pch,xlab=xlab,ylab=ylab,
                kde.col=kde.col,hist.col=hist.col,
                show.hist=show.hist,bty=bty,binwidth=binwidth,
                ncol=ncol,...)
}
#' @rdname kde
#' @export
kde.RbSr <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,log=FALSE,
                     n=512,plot=TRUE,pch=NA,xlab="age [Ma]",ylab="",
                     kde.col=rgb(1,0,1,0.6),hist.col=rgb(0,1,0,0.2),
                     show.hist=TRUE,bty='n',binwidth=NA,ncol=NA,
                     i2i=TRUE,...){
    tt <- RbSr.age(x,i2i=i2i)[,1]
    kde.default(tt,from=from,to=to,bw=bw,adaptive=adaptive,log=log,
                n=n,plot=plot,pch=pch,xlab=xlab,ylab=ylab,
                kde.col=kde.col,hist.col=hist.col,
                show.hist=show.hist,bty=bty,binwidth=binwidth,
                ncol=ncol,...)
}
#' @rdname kde
#' @export
kde.LuHf <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,log=FALSE,
                     n=512,plot=TRUE,pch=NA,xlab="age [Ma]",ylab="",
                     kde.col=rgb(1,0,1,0.6),hist.col=rgb(0,1,0,0.2),
                     show.hist=TRUE,bty='n',binwidth=NA,ncol=NA,
                     i2i=TRUE,...){
    tt <- LuHf.age(x,i2i=i2i)[,1]
    kde.default(tt,from=from,to=to,bw=bw,adaptive=adaptive,log=log,
                n=n,plot=plot,pch=pch,xlab=xlab,ylab=ylab,
                kde.col=kde.col,hist.col=hist.col,
                show.hist=show.hist,bty=bty,binwidth=binwidth,
                ncol=ncol,...)
}
#' @rdname kde
#' @export
kde.UThHe <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,log=FALSE,
                      n=512,plot=TRUE,pch=NA,xlab="age [Ma]",ylab="",
                      kde.col=rgb(1,0,1,0.6),hist.col=rgb(0,1,0,0.2),
                      show.hist=TRUE,bty='n',binwidth=NA,ncol=NA,...){
    tt <- UThHe.age(x)[,1]
    kde.default(tt,from=from,to=to,bw=bw,adaptive=adaptive,log=log,
                n=n,plot=plot,pch=pch,xlab=xlab,ylab=ylab,
                kde.col=kde.col, hist.col=hist.col,
                show.hist=show.hist,bty=bty,binwidth=binwidth,
                ncol=ncol,...)
}
#' @rdname kde
#' @export
kde.fissiontracks <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,log=FALSE,
                              n=512,plot=TRUE,pch=NA,xlab="age [Ma]",ylab="",
                              kde.col=rgb(1,0,1,0.6),hist.col=rgb(0,1,0,0.2),
                              show.hist=TRUE,bty='n',binwidth=NA,ncol=NA,...){
    tt <- fissiontrack.age(x)[,1]
    kde.default(tt,from=from,to=to,bw=bw,adaptive=adaptive,log=log,
                n=n,plot=plot,pch=pch,xlab=xlab,ylab=ylab,
                kde.col=kde.col, hist.col=hist.col,
                show.hist=show.hist,bty=bty,binwidth=binwidth,
                ncol=ncol,...)
}


# helper functions for the generic kde function
getkde <- function(x,...){ UseMethod("getkde",x) }
getkde.default <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,log=FALSE,n=512,...){
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
    out$x <- seq(from=from,to=to,length.out=n)
    if (adaptive){
        out$y <- Abramson(d,from=from,to=to,bw=bw,n=n,...)
    } else {
        out$y <- stats::density(d,bw,from=from,to=to,n=n,...)$y
    }
    if (log) out$x <- exp(out$x)
    out$y <- out$y/(sum(out$y)*(to-from)/n)
    out$x <- c(out$x[1],out$x,out$x[n])
    out$y <- c(0,out$y,0)
    out$bw <- bw
    out$ages <- x
    out
}
getkde.detritals <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,
                             log=FALSE, n=512,samebandwidth=TRUE,
                             normalise=TRUE,...){
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
                               plot=FALSE,...)
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

# Abramson
# get geometric mean pilot density
getG <- function(pdens) {
    fpos <- pdens[pdens>0]
    N <- length(fpos)
    out <- exp(sum(log(fpos))/N)
    out
}

# Abramson
# get fixed bandwidth pilot density
pilotdensity <- function(dat,bw){
    d <- stats::na.omit(dat)
    n <- length(d)
    dens <- rep(0,n)
    for (i in 1:n){
        dens[i] <- mean(stats::density(d,bw,from=(d[i]-1e-10),
                                       to=(d[i]+1e-10),n=2)$y)
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
        dens <- dens + stats::density(d[i],bw*lambda,from=from,to=to,n=n,...)$y
    }
    dens
}

plot.KDE <- function(x,pch='|',xlab="age [Ma]",ylab="",
                     kde.col=grDevices::rgb(1,0,1,0.6),show.hist=TRUE,
                     hist.col=grDevices::rgb(0,1,0,0.2),binwidth=NA,bty='n',...){
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
    graphics::points(ages,rep(graphics::par("usr")[3]/2,length(ages)),pch=pch)
    graphics::text(M,maxy,paste0("n=",length(ages)),pos=2)
}

plot.KDEs <- function(x,ncol=NA,pch=NA,xlab="age [Ma]",ylab="",
                      kde.col=grDevices::rgb(1,0,1,0.6),show.hist=TRUE,
                      hist.col=grDevices::rgb(0,1,0,0.2),binwidth=NA,bty='n',...){
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
            plot.KDE(x$kdes[[i]],pch=pch,xlab=xlab,ylab=ylab,
                     kde.col=kde.col,show.hist=show.hist,
                     hist.col=hist.col,binwidth=binwidth,
                     bty=bty,ann=FALSE,ylim=ylim,...)
            graphics::mtext(side=1,text=xlab,line=2,cex=0.8)
        } else {
            plot.KDE(x$kdes[[i]],pch=pch,xlab=xlab,ylab=ylab,
                     kde.col=kde.col,show.hist=show.hist,
                     hist.col=hist.col,binwidth=binwidth,
                     bty=bty,xaxt='n',ylim=ylim,...)
        }
        graphics::title(snames[i])
    }
    graphics::par(oldpar)
}
