#' Create a kernel density estimate
#'
#' Turns a vector of numbers into an object of class \code{KDE} using
#' a combination of the Botev (2010) bandwidth selector and the
#' Abramson (1982) adaptive kernel bandwidth modifier.
#' 
#' @param x a vector of numbers
#' @param from minimum age of the time axis. If NULL, this is set
#' automatically
#' @param to maximum age of the time axis. If NULL, this is set
#' automatically
#' @param bw the bandwidth of the KDE. If NULL, bw will be calculated
#' automatically using \code{botev()}
#' @param adaptive boolean flag controlling if the adaptive KDE
#' modifier of Abramson (1982) is used
#' @param log transform the ages to a log scale if TRUE
#' @param n horizontal resolution of the density estimate
#' @param ... optional arguments to be passed on to \code{density}
#' @return an object of class \code{KDE}, i.e. a list
#' containing the following items:
#'
#' x: horizontal plot coordinates
#'
#' y: vertical plot coordinates
#'
#' bw: the base bandwidth of the density estimate
#'
#' ages: the data values from the input to the \code{KDE} function
#' @examples
#' data(examples)
#' dens <- kde(examples$DZ[['N1']],0,3000,kernel="epanechnikov")
#' plot(dens)
#' @seealso kdes
#' @export
kde <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,log=FALSE,n=512,...){
    out <- list()
    class(out) <- "KDE"
    out$name <- deparse(substitute(x))
    out$log <- log
    if (is.na(from) | is.na(to)) {
        mM <- getmM(x,from,to,log)
        from <- mM$m
        to <- mM$M
    }
    if (log) {
        d <- log(x)
        from <- log(from)
        to <- log(to)
        bw <- bw/(stats::median(x))
    } else {
        d <- x
    }
    out$x <- seq(from=from,to=to,length.out=n)
    if (is.na(bw)){ bw <- botev(d) }
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

#' Create a list of KDEs
#'
#' Convert a list of numerical vectors into a list of objects of class
#' \code{KDE}
#' 
#' @param x a named list of vectors containing ordinal data
#' @param from minimum limit of the x-axis.
#' @param to maximum limit of the x-axis.
#' @param bw the bandwidth of the kernel density estimates. If bw =
#' NA, the bandwidth will be set automatically using \code{botev()}
#' @param samebandwidth boolean flag indicating whether the same
#' bandwidth should be used for all samples. If samebandwidth = TRUE
#' and bw = NULL, then the function will use the median bandwidth of
#' all the samples.
#' @param adaptive boolean flag switching on the adaptive bandwidth
#' modifier of Abramson (1982)
#' @param normalise boolean flag indicating whether or not the KDEs
#' should all integrate to the same value.
#' @param log boolean flag indicating whether the data should by
#' plotted on a logarithmic scale.
#' @param n horizontal resolution of the density estimates
#' @param ... optional parameters to be passed on to \code{density}
#' @return an object of class \code{KDEs}, i.e. a list containing the
#' following items:
#'
#' kdes: a named list with objects of class \code{KDE}
#'
#' from: the beginning of the common time scale
#'
#' to: the end of the common time scale
#' 
#' themax: the maximum probability density of all the KDEs
#'
#' xlabel: the x-axis label to be used by \code{plot.KDEs}
#' @examples
#' data(examples)
#' KDES <- kdes(examples$DZ,from=0,to=3000)
#' plot(KDES)
#' @seealso kde
#' @export
kdes <- function(x,from=NA,to=NA,bw=NA,samebandwidth=TRUE,
                 adaptive=TRUE,normalise=FALSE,log=FALSE,n=512,...){
    if (is.na(from) | is.na(to)) {
        mM <- getmM(unlist(x),from,to,log)
        from <- mM$m
        to <- mM$M
    }
    snames <- names(x)
    thekdes <- list()
    themax <- -1
    if (is.na(bw) & samebandwidth) bw <- commonbandwidth(x)
    for (name in snames){
        thekdes[[name]] <- kde(x[[name]],from=from,to=to,bw=bw,
                               adaptive=adaptive,log=log,n=n,...)
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
    n <- length(dat)
    dens <- rep(0,n)
    for (i in 1:n){
        dens[i] <- mean(stats::density(dat,bw,from=(dat[i]-1e-10),
                                to=(dat[i]+1e-10),n=2)$y)
    }
    dens
}

# adaptive KDE algorithm of Abramson (1982) as summarised by Jahn (2007)
Abramson <- function(dat,from,to,bw,n=512,...){
    nn <- length(dat)
    pdens <- pilotdensity(dat,bw)
    G <- getG(pdens)
    lambda <- 0
    dens <- rep(0,n)
    for (i in 1:nn){
        lambda = sqrt(G/pdens[i])
        dens <- dens + stats::density(dat[i],bw*lambda,from=from,to=to,n=n,...)$y
    }
    dens
}

#' Plot a kernel density estimate
#'
#' Plots an object of class \code{KDE}
#' 
#' @param x an object of class \code{KDE}
#' @param pch the symbol used to show the samples. May be a vector.
#'     Set \code{pch = NA} to turn them off.
#' @param xlab the label of the x-axis
#' @param ylab the label of the y-axis
#' @param kde.col the fill colour of the KDE specified as a four
#'     element vector of r, g, b, alpha values
#' @param show.hist boolean flag indicating whether a histogram should
#'     be added to the KDE
#' @param hist.col the fill colour of the histogram specified as a
#'     four element vector of r, g, b, alpha values
#' @param binwidth scalar width of the histogram bins, in Myr if
#'     \code{x$log==FALSE}, or as a fractional value if
#'     \code{x$log==TRUE}. Sturges' Rule is used if \code{binwidth==NA}
#' @param bty change to \code{"o"}, \code{"l"}, \code{"7"},
#'     \code{"c"}, \code{"u"}, or \code{"]"} if you want to draw a box
#'     around the plot
#' @param ... optional parameters to be passed on to the graphics
#'     object
#' @examples
#' data(examples)
#' dens <- kde(examples$DZ[['N1']],from=0,to=3000)
#' plot(dens)
#' @importFrom graphics axis hist par rect
#' @seealso KDE
#' @method plot KDE
#' @export
plot.KDE <- function(x,pch='|',xlab="age [Ma]",ylab="",
                     kde.col=rgb(1,0,1,0.6),show.hist=TRUE,
                     hist.col=rgb(0,1,0,0.2),binwidth=NA,bty='n',...){
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
            breaks <- exp(seq(log(m),log(M),by=binwidth))
        }
        h <- hist(log(ages),breaks=log(breaks),plot=FALSE)
    } else {
        do.log <- ''
        if (is.na(binwidth)) {
            breaks <- seq(0,M,length.out=nb+1)
        } else {
            breaks <- seq(0,M,by=binwidth)
        }
        h <- hist(ages,breaks=breaks,plot=FALSE)
    }
    nb <- length(breaks)-1
    graphics::plot(x$x,x$y,type='n',log=do.log,
                   xlab=xlab,ylab=ylab,yaxt='n',bty=bty,...)
    if (show.hist){
        rect(xleft=breaks[1:nb],xright=breaks[2:(nb+1)],
             ybottom=0,ytop=h$density,col=hist.col)
        if (par('yaxt')!='n') {
            fact <- max(h$counts)/max(h$density)
            labels <- pretty(fact*h$density)
            at <- labels/fact
            axis(2,at=at,labels=labels)
        }
    }
    graphics::polygon(x$x,x$y,col=kde.col)
    graphics::lines(x$x,x$y,col='black')
    graphics::points(ages,rep(graphics::par("usr")[3]/2,length(ages)),pch=pch)
    graphics::text(M,max(x$y),paste0("n=",length(ages)),pos=2)
}

#' Plot a list of kernel density estimates
#'
#' Plots an object of class \code{KDEs}
#' 
#' @param x an object of class \code{KDEs}
#' @param ncol scalar value indicating the number of columns over
#'     which the KDEs should be divided
#' @param pch the symbol used to show the samples. May be a vector.
#'     Set \code{pch = NA} to turn them off.
#' @param xlab the label of the x-axis
#' @param ylab the label of the y-axis
#' @param kde.col the fill colour of the KDE specified as a four
#'     element vector of r, g, b, alpha values
#' @param show.hist boolean flag indicating whether a histogram should
#'     be added to the KDE
#' @param hist.col the fill colour of the histogram specified as a
#'     four element vector of r, g, b, alpha values
#' @param binwidth scalar width of the histogram bins, in Myr if
#'     \code{x$log==FALSE}, or as a fractional value if
#'     \code{x$log==TRUE}. Sturges' Rule is used if
#'     \code{binwidth==NA}
#' @param bty change to \code{"o"}, \code{"l"}, \code{"7"},
#'     \code{"c"}, \code{"u"}, or \code{"]"} if you want to draw a box
#'     around the plot
#' @param ... optional parameters to be passed on to the graphics
#'     object
#' @importFrom graphics mtext
#' @examples
#' data(examples)
#' KDES <- kdes(examples$DZ)
#' plot(KDES)
#' @seealso KDE
#' @method plot KDEs
#' @export
plot.KDEs <- function(x,ncol=NA,pch=NA,xlab="age [Ma]",ylab="",
                      kde.col=rgb(1,0,1,0.6),show.hist=TRUE,
                      hist.col=rgb(0,1,0,0.2),binwidth=NA,bty='n',...){
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
            mtext(side=1,text=xlab,line=2,cex=0.8)
        } else {
            plot.KDE(x$kdes[[i]],pch=pch,xlab=xlab,ylab=ylab,
                     kde.col=kde.col,show.hist=show.hist,
                     hist.col=hist.col,binwidth=binwidth,
                     bty=bty,xaxt='n',ylim=ylim,...)
        }
        title(snames[i])
    }
    graphics::par(oldpar)
}

#' Plot (a) kernel density estimate(s)
#'
#' Plots geochronological datsets as kernel density estimates using a
#' combination of the Botev (2010) bandwidth selector and the Abramson
#' (1982) adaptive kernel bandwidth modifier.
#' 
#' @param x an object of class \code{UPb} or \code{detritals}
#' @param from minimum age of the time axis. If NULL, this is set
#'     automatically
#' @param to maximum age of the time axis. If NULL, this is set
#'     automatically
#' @param bw the bandwidth of the KDE. If NULL, bw will be calculated
#'     automatically using \code{botev()}
#' @param adaptive boolean flag controlling if the adaptive KDE
#'     modifier of Abramson (1982) is used
#' @param log transform the ages to a log scale if TRUE
#' @param n horizontal resolution of the density estimate
#' @param samebandwidth boolean flag indicating whether the same
#'     bandwidth should be used for all samples. If samebandwidth =
#'     TRUE and bw = NULL, then the function will use the median
#'     bandwidth of all the samples. This option is only used if
#'     \code{x} is of class \code{detritals}.
#' @param normalise boolean flag indicating whether or not the KDEs
#'     should all integrate to the same value. This option is only
#'     used if \code{x} is of class \code{detritals}.
#' @param binwidth scalar width of the histogram bins, in Myr if
#'     \code{x$log==FALSE}, or as a fractional value if
#'     \code{x$log==TRUE}. Sturges' Rule is used if
#'     \code{binwidth==NA}
#' @param ... optional arguments to be passed on to \code{plot.KDEs}
#'     (if \code{x} is of class \code{detritals} or \code{plot.KDE}
#'     otherwise
#' @examples
#' data(examples)
#' kde.plot(examples$DZ[['N2']])
#' @seealso kde kdes plot.KDE plot.KDEs
#' @export
kde.plot <- function(x,from=NA,to=NA,bw=NA,adaptive=TRUE,log=FALSE,
                     n=512,samebandwidth=TRUE,normalise=FALSE,binwidth=NA,...){
    if (methods::is(x,'detritals')){
        plot.KDEs(kdes(x,from=from,to=to,bw=bw,samebandwidth=samebandwidth,
                       adaptive=adaptive,normalise=normalise,log=log,n=n),binwidth=binwidth,...)
    } else {
        plot.KDE(kde(x,from=from,to=to,bw=bw,adaptive=adaptive,log=log,n=n),binwidth=binwidth,...)
    }
}
