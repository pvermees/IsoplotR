#' Multidimensional Scaling
#'
#' Performs classical or nonmetric Multidimensional Scaling analysis
#' @param x a dissimilarity matrix OR an object of class
#'     \code{detrital}
#' @param classical logical flag indicating whether classical
#'     (\code{TRUE}) or nonmetric (\code{FALSE}) MDS should be used
#' @param plot show the MDS configuration (if \code{shepard=FALSE}) or
#'     Shepard plot (if \code{shepard=TRUE}) on a graphical device
#' @param shepard logical flag indicating whether the graphical output
#'     should show the MDS configuration (\code{shepard=FALSE}) or a
#'     Shepard plot with the 'stress' value. This argument is only
#'     used if \code{plot=TRUE}.
#' @param nnlines if \code{TRUE}, draws nearest neighbour lines
#' @param pch plot character (see \code{?plot} for details). May be a
#'     vector.
#' @param pos a position specifier for the labels (if
#'     \code{pch!=NA}). Values of 1, 2, 3 and 4 indicate positions
#'     below, to the left of, above and to the right of the MDS
#'     coordinates, respectively.
#' @param cex.symbols a numerical value giving the amount by which
#'     plotting symbols should be magnified relative to the default
#' @param col plot colour (may be a vector)
#' @param bg background colour (may be a vector)
#' @param xlab a string with the label of the x axis
#' @param ylab a string with the label of the y axis
#' @param ... optional arguments to the generic \code{plot} function
#' @return if \code{plot=FALSE}, returns an object of class
#'     \code{MDS}, i.e. a list containing the following items:
#'     \describe{ \item{points}{a two column vector of the fitted
#'     configuration} \item{classical}{a logical flag indicating
#'     whether the MDS configuration was obtained by classical
#'     (\code{TRUE}) or nonmetric (\code{FALSE}) MDS} \item{diss}{the
#'     dissimilarity matrix used for the MDS analysis}
#'     \item{stress}{(only if \code{classical=TRUE}) the final stress
#'     achieved (in percent)} }
#' @references Vermeesch, P., 2013. Multi-sample comparison of
#'     detrital age distributions. Chemical Geology, 341, pp.140-146.
#' @examples
#' data(examples)
#' # Parameters 'xaxt' and 'yaxt' control if the axis is plotted at all.
#' mds(examples$DZ,nnlines=TRUE,cex=5,xaxt='n',yaxt='n')
#' dev.new()
#' mds(examples$DZ,shepard=TRUE)
#' @rdname mds
#' @export
mds <- function(x,...){ UseMethod("mds",x) }
#' @rdname mds
#' @export
mds.default <- function(x,classical=FALSE,plot=TRUE,shepard=FALSE,
                        nnlines=FALSE,pch=21,pos=NULL,cex.symbols=2.5,
                        col='black',bg='white',xlab="",ylab="",...){
    out <- list()
    if (classical) out$points <- stats::cmdscale(x)
    else out <- MASS::isoMDS(d=x)
    out$classical <- classical
    out$diss <- x
    class(out) <- "MDS"
    if (plot) plot.MDS(out,nnlines=nnlines,pch=pch,pos=pos,
                       cex.symbols=cex.symbols, shepard=shepard,
                       col=col,bg=bg,xlab=xlab, ylab=ylab,...)
    else return(out)
}
#' @rdname mds
#' @export
mds.detritals <- function(x,classical=FALSE,plot=TRUE,shepard=FALSE,
                        nnlines=FALSE,pch=21,pos=NULL,cex.symbols=2.5,
                        col='black',bg='white',xlab="",ylab="",...){
    d <- diss(x)
    out <- mds.default(d,classical=classical,plot=plot,
                       shepard=shepard, nnlines=nnlines,pch=pch,
                       pos=pos,cex.symbols=cex.symbols, col=col,bg=bg,
                       xlab=xlab,ylab=ylab,...)
    out
}

# x is an object of class detrital
diss <- function(x) {
    n <- length(x)
    d <- mat.or.vec(n,n)
    rownames(d) <- names(x)
    colnames(d) <- names(x)
    for (i in 1:n){
        for (j in 1:n){
            d[i,j] <- KS.diss(x[[i]],x[[j]])
    }   }
    stats::as.dist(d)
}

KS.diss <- function(x,y) {
    xx = sort(x)
    cdftmp = stats::ecdf(xx)
    cdf1 = cdftmp(xx)
    xy = sort(y)
    cdftmp = stats::ecdf(xy)
    cdfEstim = cdftmp(xy)
    cdfRef = stats::approx(xx, cdf1, xy, yleft = 0, yright = 1, ties = "mean")
    dif = cdfRef$y - cdfEstim
    dif = abs(dif)
    out = max(dif)
    return(out)
}

plot.MDS <- function(x,nnlines=FALSE,pch=21,pos=NULL,cex.symbols=2.5,
                     shepard=FALSE,col='black',bg='white',xlab="",
                     ylab="",...){
    if (shepard & !x$classical){
        shep <- MASS::Shepard(x$diss,x$points)
        graphics::plot(shep,pch=pch,col=col,bg=bg,
                       xlab='dissimilarities',
                       ylab='distances/disparities',...)
        graphics::lines(shep$x,shep$yf,type="S")
        graphics::title(paste0("Stress = ",x$stress))
    } else {
        graphics::plot(x$points,type='n',asp=1,xlab=xlab,ylab=ylab,...)
        if (nnlines) plotlines(x$points,x$diss)
        graphics::points(x$points,pch=pch,cex=cex.symbols,col=col,bg=bg)
        graphics::text(x$points,labels=labels(x$diss),pos=pos,col=col,bg=bg)
    }
}

# a function to plot the nearest neighbour lines
plotlines <- function(conf,diss) {
    # rank the samples according to their pairwise proximity
    i = t(apply(as.matrix(diss),1,function(x) order(x))[2:3,])
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
