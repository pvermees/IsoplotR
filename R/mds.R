#' @title
#' Multidimensional Scaling
#'
#' @description Performs classical or nonmetric Multidimensional
#'     Scaling analysis
#'
#' @details
#' Multidimensional Scaling (MDS) is a dimension-reducting technique
#' that takes a matrix of pairwise `dissimilarities' between objects
#' (e.g., age distributions) as input and produces a configuration of
#' two (or higher-) dimensional coordinates as output, so that the
#' Euclidean distances between these coordinates approximate the
#' dissimilarities of the input matrix. Thus, an MDS-configuration
#' serves as a `map' in which similar samples cluster closely together
#' and dissimilar samples plot far apart. In the context of detrital
#' geochronology, the dissimilarity between samples is given by the
#' statistical distance between age distributions. There are many ways
#' to define this statistical distance. \code{IsoplotR} uses the
#' Kolmogorov-Smirnov (KS) statistic due to its simplicity and the
#' fact that it behaves like a true distance in the mathematical sense
#' of the word (Vermeesch, 2013). The KS-distance is given by the
#' maximum vertical distance between two \code{\link{cad}} step
#' functions. Thus, the KS-distance takes on values between zero
#' (perfect match between two age distributions) and one (no overlap
#' between two distributions).  Calculating the KS-distance between
#' samples two at a time populates a symmetric dissimilarity matrix
#' with positive values and a zero diagonal. \code{IsoplotR}
#' implements two algorithms to convert this matrix into a
#' configuration. The first (`classical') approach uses a sequence of
#' basic matrix manipulations developed by Young and Householder
#' (1938) and Torgerson (1952) to achieve a linear fit between the
#' KS-distances and the fitted distances on the MDS configuration. The
#' second, more sophisticated (`nonmetric') approach subjects the
#' input distances to a transformation \eqn{f} prior to fitting a
#' configuration:
#' \cr\cr
#' \eqn{\delta_{i,j} = f(KS_{i,j})}
#' \cr\cr
#' where \eqn{KS_{i,j}} is the KS-distance between samples \eqn{i} and
#' \eqn{j} (for \eqn{1 \leq i \neq j \leq n}) and \eqn{\delta_{i,j}}
#' is the `disparity' (Kruskal, 1964).  Fitting an MDS
#' configuration then involves finding the disparity transformation
#' that maximises the goodness of fit (or minimises the `stress')
#' between the disparities and the fitted distances. The latter two
#' quantities can also be plotted against each other as a `Shepard
#' plot'.
#'
#' @param x a dissimilarity matrix OR an object of class
#'     \code{detrital}
#' @param method either \code{'KS'} (for the Kolmogorov-Smirnov
#'     distance) or \code{'W2'} (for the Wasserstein-2 distance).
#' @param classical logical flag indicating whether classical
#'     (\code{TRUE}) or nonmetric (\code{FALSE}) MDS should be used
#' @param plot show the MDS configuration (if \code{shepard=FALSE}) or
#'     Shepard plot (if \code{shepard=TRUE}) on a graphical device
#' @param shepard logical flag indicating whether the graphical output
#'     should show the MDS configuration (\code{shepard=FALSE}) or a
#'     Shepard plot with the 'stress' value. This argument is only
#'     used if \code{plot=TRUE}.
#' @param nnlines if \code{TRUE}, draws nearest neighbour lines
#' @param pos a position specifier for the labels (if
#'     \code{par('pch')!=NA}). Values of 1, 2, 3 and 4 indicate
#'     positions below, to the left of, above and to the right of the
#'     MDS coordinates, respectively.
#' @param col plot colour (may be a vector)
#' @param bg background colour (may be a vector)
#' @param xlab a string with the label of the x axis
#' @param ylab a string with the label of the y axis
#' @param hide vector with indices of aliquots that should be removed
#'     from the plot.
#' @param asp aspect ratio of the MDS configuration. See
#'     \code{plot.window} for further details.
#' @param ... optional arguments to the generic \code{plot} function
#' @seealso \code{\link{cad}}, \code{\link{kde}}
#' @return Returns an object of class \code{MDS}, i.e. a list
#'     containing the following items:
#'
#' \describe{
#' \item{points}{a two-column vector of the fitted configuration}
#' \item{classical}{a logical flag indicating whether the MDS
#'     configuration was obtained by classical (\code{TRUE}) or
#'     nonmetric (\code{FALSE}) MDS}
#' \item{diss}{the dissimilarity matrix used for the MDS analysis}
#' \item{stress}{(only if \code{classical=TRUE}) the final stress
#'     achieved (in percent)}
#' }
#'
#' @references
#' Kruskal, J., 1964. Multidimensional scaling by optimizing goodness
#' of fit to a nonmetric hypothesis. Psychometrika 29 (1), 1-27.
#'
#' Torgerson, W. S. Multidimensional scaling: I. Theory and
#' method. Psychometrika, 17(4): 401-419, 1952.
#'
#' Vermeesch, P., 2013. Multi-sample comparison of detrital age
#' distributions. Chemical Geology, 341, pp.140-146.
#'
#' Young, G. and Householder, A. S. Discussion of a set of points in
#' terms of their mutual distances. Psychometrika, 3(1):19-22, 1938.
#'
#' @examples
#' attach(examples)
#' mds(DZ,nnlines=TRUE,pch=21,cex=5)
#' dev.new()
#' mds(DZ,shepard=TRUE)
#' @rdname mds
#' @export
mds <- function(x,...){ UseMethod("mds",x) }
#' @rdname mds
#' @export
mds.default <- function(x,classical=FALSE,plot=TRUE,shepard=FALSE,
                        nnlines=FALSE,pos=NULL,col='black',
                        bg='white',xlab=NA,ylab=NA,asp=1,...){
    out <- list()
    if (classical) out$points <- stats::cmdscale(x)
    else out <- MASS::isoMDS(d=x)
    out$classical <- classical
    out$diss <- x
    class(out) <- "MDS"
    if (plot) plot.MDS(out,nnlines=nnlines,pos=pos,
                       shepard=shepard,col=col,bg=bg,
                       xlab=xlab,ylab=ylab,asp=asp,...)
    invisible(out)
}
#' @rdname mds
#' @export
mds.detritals <- function(x,method="KS",classical=FALSE,plot=TRUE,
                          shepard=FALSE,nnlines=FALSE,pos=NULL,col='black',
                          bg='white',xlab=NA,ylab=NA,hide=NULL,asp=1,...){
    if (is.character(hide)) hide <- which(names(x)%in%hide)
    x2plot <- clear(x,hide)
    d <- diss(x2plot,method=method)
    out <- mds.default(d,classical=classical,plot=plot,
                       shepard=shepard,nnlines=nnlines,pos=pos,
                       col=col,bg=bg,xlab=xlab,ylab=ylab,asp=asp,...)
    invisible(out)
}

#' @title Dissimilarity between detrital age distributions
#' @description Calculates the pairwise dissimilarity between detrital
#'     age distributions, using either the Wasserstein-2 or
#'     Kolmogorov-Smirnov distance.
#' @details The Kolmogorov-Smirnov statistic is the maximum vertical
#'     difference between two empirical cumulative distribution
#'     functions. The Wasserstein distance is a function of the area
#'     between them. Both dissimilarity measures are useful for
#'     multidimensional scaling.
#' @param x an object of class \code{detrital} OR a vector of numbers
#' @param y a vector of numbers
#' @param method either \code{'KS'} (for Kolmogorov-Smirnov distance),
#'     or \code{'W2'} (for Wasserstein-2 distance).
#' @author Written by Pieter Vermeesch, using modified code from
#'     Mathieu Vrac's \code{CDFt} package (\code{KolmogorovSmirnov}
#'     function), and Dominic Schuhmacher's \code{transport} package
#'     (\code{transport1d} function).
#' @seealso \code{\link{mds}}
#' @param ... extra arguments (not used)
#' @return an object of class \code{dist}.
#' @examples
#' d <- diss(examples$DZ,method='KS')
#' mds(d)
#' @rdname diss
#' @export
diss <- function(x,...){ UseMethod("diss",x) }
#' @rdname diss
#' @export
diss.default <- function(x,y,method='KS',...){
    if (identical(method,'W2')){
        out <- Wasserstein.diss(x,y)
    } else {
        out <- KS.diss(x,y)
    }
    out
}
#' @rdname diss
#' @export
diss.detritals <- function(x,method='W2',...) {
    n <- length(x)
    d <- mat.or.vec(n,n)
    rownames(d) <- names(x)
    colnames(d) <- names(x)
    for (i in 1:n){
        for (j in 1:n){
            d[i,j] <- diss(x[[i]],x[[j]],method=method,...)
    }   }
    stats::as.dist(d)
}
KS.diss <- function(x,y) {
    xx <- sort(x)
    cdftmp <- stats::ecdf(xx)
    cdf1 <- cdftmp(xx)
    xy <- sort(y)
    cdftmp <- stats::ecdf(xy)
    cdfEstim <- cdftmp(xy)
    cdfRef <- stats::approx(xx,cdf1,xy,yleft=0,yright=1,ties="mean")
    dif <- abs(cdfRef$y - cdfEstim)
    max(dif)
}
# modified after the wasserstein1d function of the transport package
Wasserstein.diss <- function(x,y,p=2) {
    m <- length(x)
    n <- length(y)
    stopifnot(m>0 && n>0)
    if (m == n) {
        out <- mean(abs(sort(y)-sort(x))^p)^(1/p)
    } else {
        wx <- rep(1,m)
        wy <- rep(1,n)
        ordx <- order(x)
        ordy <- order(y)
        x <- x[ordx]
        y <- y[ordy]
        wx <- wx[ordx]
        wy <- wy[ordy]
        ux <- (wx/sum(wx))[-m]
        uy <- (wy/sum(wy))[-n]
        cux <- c(cumsum(ux))
        cuy <- c(cumsum(uy))
        xrep <- graphics::hist(cuy,breaks=c(-Inf,cux,Inf),plot=FALSE)$counts+1
        yrep <- graphics::hist(cux,breaks=c(-Inf,cuy,Inf),plot=FALSE)$counts+1
        xx <- rep(x,times=xrep)
        yy <- rep(y,times=yrep)
        uu <- sort(c(cux,cuy))
        uu0 <- c(0,uu)
        uu1 <- c(uu,1)
        out <- sum((uu1-uu0)*abs(yy-xx)^p)^(1/p)
    }
    out
}

plot.MDS <- function(x,nnlines=FALSE,pos=NULL,shepard=FALSE,
                     col='black',bg='white',xlab=NA,ylab=NA,asp=1,pch,...){
    if (shepard & !x$classical){
        if (missing(pch)) pch <- 21
        if (is.na(xlab)) xlab <- 'dissimilarities'
        if (is.na(ylab)) ylab <- 'distances/disparities'
        shep <- MASS::Shepard(x$diss,x$points)
        graphics::plot(x=shep,col=col,bg=bg,pch=pch,
                       xlab=xlab,ylab=ylab,...)
        graphics::lines(shep$x,shep$yf,type="S")
        graphics::title(paste0("Stress = ",x$stress))
    } else {
        if (missing(pch)) pch <- NA
        if (is.na(xlab)) xlab <- 'Dim 1'
        if (is.na(ylab)) ylab <- 'Dim 2'
        graphics::plot(x=x$points,type='n',asp=asp,
                       xlab=xlab,ylab=ylab,...)
        if (nnlines) plotlines(x$points,x$diss)
        if (is.na(pch)) {
            graphics::points(x$points,pch=pch,...)
            graphics::text(x=x$points,labels=labels(x$diss),
                           col=col,bg=bg,pos=pos)
        } else {
            graphics::points(x=x$points,pch=pch,col=col,bg=bg,...)
            graphics::text(x$points,labels=labels(x$diss),pos=pos)
        }
    }
}

# a function to plot the nearest neighbour lines
plotlines <- function(conf,diss) {
    # rank the samples according to their pairwise proximity
    i = t(apply(as.matrix(diss),1,function(x) order(x))[2:3,,drop=FALSE])
    # coordinates for the lines
    x1 = as.vector(conf[i[,1],1]) # calculate (x,y)-coordinates ...
    y1 = as.vector(conf[i[,1],2]) # ... of nearest neighbours
    x2 = as.vector(conf[i[,2],1]) # calculate (x,y)-coordinates ...
    y2 = as.vector(conf[i[,2],2]) # ... of second nearest neighbours
    for (j in 1:nrow(conf)) {
        graphics::lines(c(conf[j,1],x1[j]),c(conf[j,2],y1[j]),lty=1) # solid line
        graphics::lines(c(conf[j,1],x2[j]),c(conf[j,2],y2[j]),lty=2) # dashed line
    }
}
