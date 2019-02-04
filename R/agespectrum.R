#' Plot a (40Ar/39Ar) release spectrum
#'
#' Produces a plot of boxes whose widths correspond to the cumulative
#' amount of \eqn{^{39}}Ar (or any other variable), and whose
#' heights express the analytical uncertainties.  Only propagates the
#' analytical uncertainty associated with decay constants and
#' J-factors \emph{after} computing the plateau composition.
#'
#' @details
#' \code{IsoplotR} defines the `plateau age' as the weighted mean age
#' of the longest sequence (in terms of cumulative \eqn{^{39}}Ar
#' content) of consecutive heating steps that pass the modified
#' Chauvenet criterion (see \code{\link{weightedmean}}).  Note that
#' this definition is different (and simpler) than the one used by
#' \code{Isoplot} (Ludwig, 2003). However, it is important to mention
#' that all definitions of an age plateau are heuristic by nature and
#' should not be used for quantitative inference.
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
#' @param alpha the confidence level of the error bars/boxes and
#'     confidence intervals.
#' @param plateau logical flag indicating whether a plateau age should
#'     be calculated. If \code{plateau=TRUE}, the function will
#'     compute the weighted mean of the largest succession of steps
#'     that pass the Chi-square test for age homogeneity.  If
#'     \code{TRUE}, returns a list with plateau parameters.
#' @param levels a vector with additional values to be displayed as
#'     different background colours of the plot symbols.
#' @param plateau.col a vector of two fill colours of the rectangles
#'     used to mark the steps belonging to the age plateau.  If
#'     \code{levels=NA}, then only the first colour is used. If
#'     \code{levels} is a vector of numbers, then \code{bg} is used to
#'     construct a colour ramp.
#' @param non.plateau.col if \code{plateau=TRUE}, the steps that do
#'     NOT belong to the plateau are given a different colour.
#' @param clabel label of the colour legend
#' @param sigdig the number of significant digits of the numerical
#'     values reported in the title of the graphical output.
#' @param line.col colour of the average age line
#' @param lwd width of the average age line
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param random.effects if \code{TRUE}, computes the weighted mean
#'     using a random effects model with two parameters: the mean and
#'     the dispersion. This is akin to a `model-3' isochron
#'     regression.
#' 
#'     if \code{FALSE}, attributes any excess dispersion to an
#'     underestimation of the analytical uncertainties. This akin to a
#'     `model-1' isochron regression.
#' @param hide vector with indices of aliquots that should be removed
#'     from the plot.
#' @param ... optional parameters to the generic \code{plot} function
#'
#' @return If \code{plateau=TRUE}, returns a list with the following
#'     items:
#'
#' \describe{
#' \item{mean}{a 3-element vector with:
#'
#' \code{x}: the plateau mean
#'
#' \code{s[x]}: the standard error of \code{x}
#'
#' \code{ci[x]}: the width of a 100(1-\eqn{\alpha})\% confidence interval of
#' \code{t} }
#'
#' \item{disp}{a 3-element vector with:
#'
#' \code{w}: the overdispersion, i.e. the standard deviation of the
#' Normal distribution that is assumed to describe the true ages.
#'
#' \code{ll}: the width of the lower half of a 100(1-\eqn{\alpha})\%
#' confidence interval for the overdispersion
#'
#' \code{ul}: the width of the upper half of a 100(1-\eqn{\alpha})\%
#' confidence interval for the overdispersion}
#'
#' \item{df}{the degrees of freedom for the weighted mean plateau fit}
#'
#' \item{mswd}{the mean square of the weighted deviates of the plateau}
#'
#' \item{p.value}{the p-value of a Chi-square test with \eqn{df=n-2}
#' degrees of freedom, where \eqn{n} is the number of steps in the
#' plateau and 2 degrees of freedom have been removed to estimate the
#' mean and the dispersion.}
#'
#' \item{fract}{the fraction of \eqn{^{39}}Ar contained in the
#' plateau}
#'
#' \item{plotpar}{plot parameters for the weighted mean (see
#' \code{\link{weightedmean}})}
#'
#' \item{i}{indices of the steps that are retained for the plateau age
#' calculation}
#'
#' }
#'
#' @seealso \code{\link{weightedmean}}
#' @rdname agespectrum
#' @export
agespectrum <- function(x,...){ UseMethod("agespectrum",x) }
#' @importFrom grDevices rgb
#' @rdname agespectrum
#' @export
agespectrum.default <- function(x,alpha=0.05,plateau=TRUE,
                                random.effects=TRUE,levels=NA,clabel="",
                                plateau.col=c("#00FF0080","#FF000080"),
                                non.plateau.col="#00FFFF80",
                                sigdig=2,line.col='red',lwd=2,
                                xlab='cumulative fraction',
                                ylab='age [Ma]',hide=NULL,...){
    XY <- plot.spectrum.axes(x=x,alpha=alpha,xlab=xlab,
                             ylab=ylab,hide=hide,levels=levels,
                             plateau.col=plateau.col,clabel=clabel,...)
    pc <- get.plateau.colours(x=x,levels=levels,plateau=plateau,
                              hide=hide,plateau.col=plateau.col,
                              non.plateau.col=non.plateau.col,
                              random.effects=random.effects,alpha=alpha)
    if (plateau){
        plot.plateau(fit=pc$plat,line.col=line.col,lwd=lwd)
        graphics::title(plateau.title(pc$plat,sigdig=sigdig,Ar=FALSE))
    }
    plot.spectrum(XY=XY,col=pc$col)
    if (plateau) return(invisible(pc$plat))
}
#' @param i2i `isochron to intercept':
#'     calculates the initial (aka `inherited',
#'     `excess', or `common') \eqn{^{40}}Ar/\eqn{^{36}}Ar ratio from
#'     an isochron fit. Setting \code{i2i} to \code{FALSE} uses the
#'     default values stored in \code{settings('iratio',...)}
#' @param exterr propagate the external (decay constant and
#'     calibration factor) uncertainties?
#' @examples
#' data(examples)
#' agespectrum(examples$ArAr,ylim=c(0,80))
#' @rdname agespectrum
#' @export
agespectrum.ArAr <- function(x,alpha=0.05,plateau=TRUE,
                             random.effects=TRUE,levels=NA,clabel="",
                             plateau.col=c("#00FF0080","#FF000080"),
                             non.plateau.col="#00FFFF80",sigdig=2,
                             exterr=TRUE,line.col='red',lwd=2,
                             i2i=FALSE,hide=NULL,...){
    x <- clear(x,hide)
    tt <- ArAr.age(x,exterr=FALSE,i2i=i2i)
    X <- cbind(x$x[,'Ar39'],tt)
    x.lab <- expression(paste("cumulative ",""^"39","Ar fraction"))
    y.lab='age [Ma]'
    XY <- plot.spectrum.axes(x=X,alpha=alpha,xlab=x.lab,
                             ylab=y.lab,hide=hide,levels=levels,
                             plateau.col=plateau.col,clabel=clabel,...)
    pc <- get.plateau.colours(x=X,levels=levels,plateau=plateau,
                              hide=hide,plateau.col=plateau.col,
                              non.plateau.col=non.plateau.col,
                              random.effects=random.effects,alpha=alpha)
    if (plateau){
        if (exterr) pc$plat <- add.exterr.to.wtdmean(x,pc$plat)
        plot.plateau(fit=pc$plat,line.col=line.col,lwd=lwd)
        graphics::title(plateau.title(pc$plat,sigdig=sigdig,
                                      Ar=TRUE,units='Ma'))
    }
    plot.spectrum(XY=XY,col=pc$col)
    if (plateau) return(invisible(pc$plat))
}

plot.spectrum.axes <- function(x,alpha=0.05,xlab='cumulative fraction',
                               ylab='age [Ma]',hide=NULL,levels=NA,
                               plateau.col=c("#00FF0080","#FF000080"),
                               clabel="",...){
    ns <- nrow(x)
    x <- clear(x[,1:3],hide)
    valid <- !is.na(rowSums(x))
    X <- c(0,cumsum(x[valid,1])/sum(x[valid,1]))
    Y <- x[valid,2]
    sY <- x[valid,3]
    fact <- stats::qnorm(1-alpha/2)
    Yl <- Y-fact*sY
    Yu <- Y+fact*sY
    minY <- min(Yl,na.rm=TRUE)
    maxY <- max(Yu,na.rm=TRUE)
    graphics::plot(c(0,1),c(minY,maxY),type='n',xlab=xlab,ylab=ylab,...)
    colourbar(z=levels,col=plateau.col,clabel=clabel)
    list(X=X,Yl=Yl,Yu=Yu,ylim=c(minY,maxY))
}
get.plateau.colours <- function(x,levels=NA,plateau=TRUE,hide=NULL,
                                plateau.col=c("#00FF0080","#FF000080"),
                                non.plateau.col="#00FFFF80",
                                random.effects=TRUE,alpha=0.05){
    ns <- nrow(x)
    if (!all(is.na(levels))) levels <- clear(levels,hide)
    if (plateau){
        plat <- get.plateau(x,alpha=alpha,random.effects=random.effects)
        plat$valid <- NULL
        colour <- rep(non.plateau.col,ns)
        np <- length(plat$i)
        levels <- levels[plat$i]
        cols <- set.ellipse.colours(ns=np,levels=levels,col=plateau.col)
        colour[plat$i] <- cols
        plat$n <- ns
    } else {
        plat <- NA
        colour <- set.ellipse.colours(ns=ns,levels=levels,
                                      hide=hide,col=plateau.col)
    }
    list(col=colour,plat=plat)
}
plot.spectrum <- function(XY,col){
    ns <- length(XY$X)
    for (i in 1:ns){
        graphics::rect(XY$X[i],XY$Yl[i],XY$X[i+1],XY$Yu[i],col=col[i])
        if (i<ns) graphics::lines(rep(XY$X[i+1],2),c(XY$Yl[i],XY$Yu[i+1]))
    }
}
plot.plateau <- function(fit,line.col='red',lwd=2){
    ci.exterr <- fit$plotpar$ci.exterr
    if (!all(is.na(ci.exterr))){
        ci.exterr$x <- c(0,1,1,0)
        graphics::polygon(ci.exterr,col='gray90',border=NA)
    }
    ci <- fit$plotpar$ci
    ci$x <- c(0,1,1,0)
    graphics::polygon(ci,col='gray75',border=NA)
    graphics::lines(c(0,1),rep(fit$mean[1],2),col=line.col,lwd=lwd)
}
plateau.title <- function(fit,sigdig=2,Ar=TRUE,units='',...){
    rounded.mean <- roundit(fit$mean['x'],
                            fit$mean[c('s[x]','ci[x]')],
                            sigdig=sigdig)
    line1 <- substitute('mean ='~a%+-%b~'|'~c~u~'(n='*n/N*')',
                        list(a=rounded.mean[1],
                             b=rounded.mean[2],
                             c=rounded.mean[3],
                             u=units,
                             n=length(fit$i),
                             N=fit$n))
    a <- signif(100*fit$fract,sigdig)
    if (Ar) line2 <- bquote(paste("Includes ",.(a),"% of the ",""^"39","Ar"))
    else line2 <- bquote(paste("Includes ",.(a),"% of the spectrum"))
    mymtext(line1,line=1,...)
    mymtext(line2,line=0,...)
}
# x is a three column vector with Ar39
# cumulative fractions, ages and uncertainties
get.plateau <- function(x,alpha=0.05,random.effects=TRUE){
    X <- x[,1]/sum(x[,1],na.rm=TRUE)
    YsY <- subset(x,select=c(2,3))
    ns <- length(X)
    out <- list()
    out$mean <- YsY[1,1:2]
    out$mswd <- 1
    out$p.value <- 1
    out$fract <- 0
    for (i in 1:(ns-1)){ # at least two steps in a plateau
        for (j in ns:(i+1)){
            fract <- sum(X[i:j],na.rm=TRUE)
            Y <- YsY[i:j,1]
            sY <- YsY[i:j,2]
            valid <- chauvenet(Y,sY,valid=rep(TRUE,j-i+1),
                               random.effects=random.effects)
            if (all(valid) & (fract > out$fract)){
                out <- weightedmean(YsY[i:j,],random.effects=random.effects,
                                    plot=FALSE,detect.outliers=FALSE,
                                    alpha=alpha)
                out$i <- i:j
                out$fract <- fract
                break
            }
        }
    }
    out
}
