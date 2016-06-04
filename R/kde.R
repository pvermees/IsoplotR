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
#' data(DZ)
#' dens <- kde(DZ[['N1']],0,3000,kernel="epanechnikov")
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
#' @param pch (optional) symbol to be used to mark the sample points along the x-axis
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
#' pch: the plot symbol to be used by \code{plot.KDEs}
#'
#' xlabel: the x-axis label to be used by \code{plot.KDEs}
#' @examples
#' data(DZ)
#' KDES <- kdes(DZ,0,3000,pch=NA)
#' summaryplot(KDES,ncol=3)
#' @seealso kde
#' @export
kdes <- function(x,from=NA,to=NA,bw=NA,samebandwidth=TRUE,
                 adaptive=TRUE,pch=NA,normalise=FALSE,
                 log=FALSE,n=512,...){
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
    out$pch <- pch
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
#' @param ... optional parameters to be passed on to the graphics
#'     object
#' @examples
#' data(Namib)
#' samp <- Namib$DZ$x[['N1']]
#' dens <- KDE(samp,from=0,to=3000)
#' plot(dens)
#' @seealso KDE
#' @method plot KDE
#' @export
plot.KDE <- function(x,pch='|',xlab="age [Ma]",ylab="",kde.col=rgb(1,0,1,0.6),
                     show.hist=TRUE,hist.col=rgb(0,1,0,0.2),
                     binwidth=NA,bty='n',...){
    if (is.na(binwidth)) breaks <- "Sturges"
    if (x$log){
        do.log <- 'x'
        if (!is.na(binwidth))
            breaks <- exp(seq(log(x$x[1]),log(tail(x$x,n=1)),by=binwidth))
        h <- hist(log(x$ages),breaks=breaks,plot=FALSE)
        width <- diff(exp(h$breaks))
    } else {
        do.log <- ''
        if (!is.na(binwidth)) breaks <- seq(0,tail(x$x,1),by=binwidth)
        h <- hist(x$ages,breaks=breaks,plot=FALSE)
        width <- diff(h$breaks)
    }
    graphics::plot(x$x,x$y,type='n',log=do.log,xlab=xlab,ylab=ylab,yaxt='n',bty=bty,...)
    if (show.hist){
        barplot(height=h$density,width=width,
                space=0,add=TRUE,col=hist.col,yaxt='n',...)
        if (par('yaxt')!='n') {
            fact <- max(h$counts)/max(h$density)
            labels <- pretty(fact*h$density)
            at <- labels/fact
            axis(2,at=at,labels=labels)
        }
    }
    graphics::polygon(x$x,x$y,col=kde.col)
    graphics::lines(x$x,x$y,col='black')
    graphics::points(x$ages,rep(graphics::par("usr")[3]/2,length(x$ages)),pch=pch)
    graphics::text(utils::tail(x$x,n=1),.9*max(x$y),paste0("n=",length(x$ages)),pos=2)
}

plot.KDEs <- function(x,ncol=1,...){
    oldpar <- graphics::par(no.readonly=T)
    snames <- names(x$kdes)
    ns <- length(snames)
    w <- rep(1,ncol) # column widths
    nppc <- ceiling(ns/ncol)
    np <- nppc*ncol # number of subpanels
    graphics::layout(matrix(1:np,nppc,length(w)),w,rep(1,nppc))
    si <- ceiling(seq(from=0,to=ns,length.out=ncol+1)) # sample index
    graphics::par(xpd=TRUE,mar=rep(1,4),oma=c(3,1,1,1))#, mfcol=c(nppc,nd))
    for (i in 1:ns){
        if ((i%%nppc)==0 | (i==ns)) plot.KDE(x$kdes[[i]],...)
        else plot.KDE(x$kdes[[i]],xaxt='n',...)
        title(snames[i])
    }
    graphics::par(oldpar)
}
