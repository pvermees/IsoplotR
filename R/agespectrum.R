#' Plot a (\eqn{^{40}}Ar/\eqn{^{39}}Ar) release spectrum
#'
#' Produces a plot of boxes whose widths correspond to the cumulative
#' amount of \eqn{^{39}}Ar (or any other volume proxy), and whose
#' heights express the analytical uncertainties.
#'
#' @param x a three column matrix whose first column gives the amount
#'     of \eqn{^{39}}Ar in each aliquot, and whose second and third
#'     columns give the age and its uncertainty.
#'
#' OR
#'
#' an object of class \code{ArAr} with \code{format=2}
#' 
#' @param alpha the confidence limits of the error bars/boxes.
#' @param plateau Boolean flag indicating whether a plateau age should
#'     be calculated if \code{plateau=TRUE}, the function will compute
#'     the weighted mean of the largest succession of steps that yield
#'     values passing the Chi-square test for age homogeneity.
#' @param plateau.col the fill colour of the rectangles used to mark
#'     the steps belonging to the age plateau.
#' @param non.plateau.col if \code{plateau=TRUE}, the steps that do
#'     NOT belong to the plateau are given a different colour.
#' @param sigdig the number of significant digits of the numerical
#'     values reported in the title of the graphical output
#'     (only used if \code{plateau=FALSE}).
#' @param line.col colour of the isochron line
#' @param lwd line width
#' @rdname agespectrum
#' @export
agespectrum <- function(x,...){ UseMethod("agespectrum",x) }
#' @rdname agespectrum
#' @export
agespectrum.default <- function(x,alpha=0.05,plateau=TRUE,
                                plateau.col=rgb(0,1,0,0.5),
                                non.plateau.col=rgb(0,1,1,0.5),
                                sigdig=2,line.col='red',lwd=2,...){
    ns <- nrow(x)
    X <- c(0,cumsum(x[,1])/sum(x[,1]))
    Y <- x[,2]
    sY <- x[,3]
    fact <- stats::qnorm(1-alpha/2)
    maxY <- max(Y+fact*sY)
    minY <- min(Y-fact*sY)
    plot(c(0,1),c(minY,maxY),type='n',...)
    if (plateau) {
        plat <- plateau(x,alpha=alpha)
        colour <- rep(non.plateau.col,ns)
        colour[plat$i] <- plateau.col
        title(plateau.title(plat,sigdig=sigdig))
        lines(c(0,1),rep(plat$mean[1],2),col=line.col,lwd=lwd)
    } else {
        colour <- rep(plateau.col,ns)
    }
    for (i in 1:ns){
        rect(X[i],Y[i]-fact*sY[i],X[i+1],Y[i]+fact*sY[i],col=colour[i])
        if (i<ns) lines(rep(X[i+1],2),c(Y[i]-fact*sY[i],Y[i+1]+fact*sY[i+1]))
    }
}
#' @param dcu propagate the decay constant uncertainties?
#' @examples
#' data(examples)
#' agespectrum(examples$ArAr)
#' @rdname agespectrum
#' @export
agespectrum.ArAr <- function(x,alpha=0.05,plateau=TRUE,
                             plateau.col=rgb(0,1,0,0.5),
                             non.plateau.col=rgb(0,1,1,0.5),sigdig=2,
                             dcu=TRUE,line.col='red',lwd=2,...){
    tt <- age(x,dcu=FALSE)
    if (x$format==2){
        X <- cbind(x$x[,'Ar39'],tt)
    } else {
        X <- cbind(seq(0,1,length.out=nrow(x$x)),tt)
    }
    x.lab <- expression(paste("cumulative ",""^"39","Ar fraction"))
    agespectrum.default(X,alpha=alpha,xlab=x.lab,ylab='age [Ma]',
                        plateau=TRUE,sigdig=sigdig,line.col=line.col,
                        lwd=lwd,...)
}

# x is a three column vector with Ar39 cumulative fractions, ages and uncertainties
plateau <- function(x,alpha=0.05){
    X <- x[,1]/sum(x[,1])
    YsY <- x[,c(2,3)]
    ns <- length(X)
    out <- list()
    out$mean <- c(0,0)
    out$mswd <- 0
    out$p.value <- 0
    out$fract <- 0
    for (i in 1:(ns-1)){
        for (j in (i+1):ns){
            fract <- sum(X[i:j])
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
                        list(a=rounded.mean$x, b=rounded.mean$err))
    line2 <- substitute('MSWD ='~a~', p('~chi^2*')='~b,
                        list(a=signif(fit$mswd,2),
                             b=signif(fit$p.value,2)))
    a <- signif(100*fit$fract,sigdig)
    line3 <- bquote(paste("Includes ",.(a),"% of the",""^"39","Ar"))
    graphics::mtext(line1,line=2)
    graphics::mtext(line2,line=1)
    graphics::mtext(line3,line=0)
}
