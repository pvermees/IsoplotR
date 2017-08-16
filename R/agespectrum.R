#' Plot a (40Ar/39Ar) release spectrum
#'
#' Produces a plot of boxes whose widths correspond to the cumulative
#' amount of \eqn{^{39}}Ar (or any other volume proxy), and whose
#' heights express the analytical uncertainties.  Only propagates the
#' analytical uncertainty associated with decay constants and
#' J-factors after computing the plateau composition.
#'
#' @param x
#' a three-column matrix whose first column gives the amount
#'     of \eqn{^{39}}Ar in each aliquot, and whose second and third
#'     columns give the age and its uncertainty.
#' 
#' OR
#' 
#' an object of class \code{ArAr}
#' 
#' @param alpha the confidence limits of the error bars/boxes.
#' @param plateau logical flag indicating whether a plateau age should
#'     be calculated. If \code{plateau=TRUE}, the function will
#'     compute the weighted mean of the largest succession of steps
#'     that yield values passing the Chi-square test for age
#'     homogeneity.
#' @param plateau.col the fill colour of the rectangles used to mark
#'     the steps belonging to the age plateau.
#' @param non.plateau.col if \code{plateau=TRUE}, the steps that do
#'     NOT belong to the plateau are given a different colour.
#' @param sigdig the number of significant digits of the numerical
#'     values reported in the title of the graphical output (only used
#'     if \code{plateau=TRUE}).
#' @param line.col colour of the isochron line
#' @param lwd line width
#' @param title add a title to the plot? If \code{FALSE}, returns a
#'     list with plateau parameters.
#' @param ... optional parameters to the generic \code{plot} function
#' @return if \code{title=FALSE}, returns a list with the following
#'     items:
#'
#' \describe{
#' \item{mean}{a 2-element vector with the plateau mean and standard error}
#'
#' \item{mswd}{the mean square of the weighted deviates of the plateau}
#'
#' \item{p.value}{the p-value of a Chi-square test with \eqn{n-1}
#' degrees of freedom, where \eqn{n} is the number of steps in the
#' plateau.}
#'
#' \item{fract}{the fraction of \eqn{^{39}}Ar contained in the
#' plateau} }
#' @rdname agespectrum
#' @export
agespectrum <- function(x,...){ UseMethod("agespectrum",x) }
#' @importFrom grDevices rgb
#' @rdname agespectrum
#' @export
agespectrum.default <- function(x,alpha=0.05,plateau=TRUE,
                                plateau.col=rgb(0,1,0,0.5),
                                non.plateau.col=rgb(0,1,1,0.5),
                                sigdig=2,line.col='red',lwd=2,
                                title=TRUE,...){
    ns <- nrow(x)
    valid <- !is.na(rowSums(x))
    X <- c(0,cumsum(x[valid,1])/sum(x[valid,1]))
    Y <- x[valid,2]
    sY <- x[valid,3]
    fact <- stats::qnorm(1-alpha/2)
    maxY <- max(Y+fact*sY,na.rm=TRUE)
    minY <- min(Y-fact*sY,na.rm=TRUE)
    graphics::plot(c(0,1),c(minY,maxY),type='n',...)
    plat <- plateau(x,alpha=alpha)
    if (plateau) {
        colour <- rep(non.plateau.col,ns)
        colour[plat$i] <- plateau.col
        graphics::lines(c(0,1),rep(plat$mean[1],2),col=line.col,lwd=lwd)
    } else {
        colour <- rep(plateau.col,ns)
    }
    for (i in 1:ns){
        graphics::rect(X[i],Y[i]-fact*sY[i],X[i+1],Y[i]+fact*sY[i],col=colour[i])
        if (i<ns) graphics::lines(rep(X[i+1],2),c(Y[i]-fact*sY[i],Y[i+1]+fact*sY[i+1]))
    }
    if (plateau & title) graphics::title(plateau.title(plat,sigdig=sigdig))
    else return(plat)
}
#' @param i2i `isochron to intercept': calculates the initial (aka `inherited',
#'     `excess', or `common') \eqn{^{40}}Ar/\eqn{^{36}}Ar ratio from an
#'     isochron fit. Setting \code{i2i} to \code{FALSE} uses the
#'     default values stored in \code{settings('iratio',...)}
#' @param exterr propagate the external (decay constant and
#'     calibration factor) uncertainties?
#' @examples
#' data(examples)
#' agespectrum(examples$ArAr,ylim=c(0,80))
#' @rdname agespectrum
#' @export
agespectrum.ArAr <- function(x,alpha=0.05,plateau=TRUE,
                             plateau.col=rgb(0,1,0,0.5),
                             non.plateau.col=rgb(0,1,1,0.5),sigdig=2,
                             exterr=TRUE,line.col='red',lwd=2,
                             i2i=FALSE,...){
    tt <- ArAr.age(x,jcu=FALSE,exterr=FALSE,i2i=i2i)
    X <- cbind(x$x[,'Ar39'],tt)
    x.lab <- expression(paste("cumulative ",""^"39","Ar fraction"))
    plat <- agespectrum.default(X,alpha=alpha,xlab=x.lab,ylab='age [Ma]',
                                plateau=plateau,sigdig=sigdig,line.col=line.col,
                                lwd=lwd,title=FALSE,...)
    # calculate the weighted mean Ar40Ar39 ratio from the weighted mean age
    R <- get.ArAr.ratio(plat$mean[1],plat$mean[2],x$J[1],0,exterr=FALSE)
    # recalculate the weighted mean age, this time
    # taking into account decay and J uncertainties
    plat$mean <- get.ArAr.age(R[1],R[2],x$J[1],x$J[2],exterr=exterr)
    if (plateau){
        graphics::title(plateau.title(plat,sigdig=sigdig))
    }
}

# x is a three column vector with Ar39 cumulative fractions, ages and uncertainties
plateau <- function(x,alpha=0.05){
    X <- x[,1]/sum(x[,1],na.rm=TRUE)
    YsY <- x[,c(2,3)]
    ns <- length(X)
    out <- list()
    out$mean <- c(0,0)
    out$mswd <- 0
    out$p.value <- 0
    out$fract <- 0
    for (i in 1:(ns-1)){
        for (j in (i+1):ns){
            fract <- sum(X[i:j],na.rm=TRUE)
            avg <- weightedmean(YsY[i:j,],plot=FALSE,
                                detect.outliers=FALSE)
            if (avg$p.value < alpha) {
                break
            } else if (fract > out$fract) {
                out$i <- i:j
                out$mean <- avg$mean
                out$mswd <- avg$mswd
                out$p.value <- avg$p.value
                out$fract <- fract
            }
        }
    }
    out
}

plateau.title <- function(fit,sigdig=2){
    rounded.mean <- roundit(fit$mean[1],fit$mean[2],sigdig=sigdig)
    line1 <- substitute('mean ='~a%+-%b~' (1'~sigma~')',
                        list(a=rounded.mean[1], b=rounded.mean[2]))
    line2 <- substitute('MSWD ='~a~', p('~chi^2*')='~b,
                        list(a=signif(fit$mswd,sigdig),
                             b=signif(fit$p.value,sigdig)))
    a <- signif(100*fit$fract,sigdig)
    line3 <- bquote(paste("Includes ",.(a),"% of the",""^"39","Ar"))
    graphics::mtext(line1,line=2)
    graphics::mtext(line2,line=1)
    graphics::mtext(line3,line=0)
}
