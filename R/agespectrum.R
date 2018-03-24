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
#' @param plateau.col the fill colour of the rectangles used to mark
#'     the steps belonging to the age plateau.
#' @param non.plateau.col if \code{plateau=TRUE}, the steps that do
#'     NOT belong to the plateau are given a different colour.
#' @param sigdig the number of significant digits of the numerical
#'     values reported in the title of the graphical output (only used
#'     if \code{plateau=TRUE}).
#' @param line.col colour of the average age line
#' @param lwd width of the average age line
#' @param title add a title to the plot?
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param show.ci show a 100(1-\eqn{\alpha})\% confidence interval for
#'     the plateau age as a grey band
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
#' \code{s[x]}: the estimated standard deviation of \code{x}
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
#' \code{\link{weightedmean}}), which are not used in the age
#' spectrum}
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
                                plateau.col=rgb(0,1,0,0.5),
                                non.plateau.col=rgb(0,1,1,0.5),
                                sigdig=2,line.col='red', lwd=2,
                                title=TRUE,show.ci=TRUE,
                                xlab='cumulative fraction',
                                ylab='age [Ma]',...){
    ns <- nrow(x)
    valid <- !is.na(rowSums(x))
    X <- c(0,cumsum(x[valid,1])/sum(x[valid,1]))
    Y <- x[valid,2]
    sY <- x[valid,3]
    fact <- stats::qnorm(1-alpha/2)
    maxY <- max(Y+fact*sY,na.rm=TRUE)
    minY <- min(Y-fact*sY,na.rm=TRUE)
    graphics::plot(c(0,1),c(minY,maxY),type='n',xlab=xlab,ylab=ylab,...)
    plat <- plateau(x,alpha=alpha)
    plat$valid <- NULL
    if (plateau) {
        colour <- rep(non.plateau.col,ns)
        colour[plat$i] <- plateau.col
        if (show.ci){
            ci <- plat$plotpar$rect
            ci$x <- c(0,1,1,0)
            graphics::polygon(ci,col='gray80',border=NA)
        }
        graphics::lines(c(0,1),rep(plat$mean[1],2),col=line.col,lwd=lwd)
    } else {
        colour <- rep(plateau.col,ns)
    }
    for (i in 1:ns){
        graphics::rect(X[i],Y[i]-fact*sY[i],
                       X[i+1],Y[i]+fact*sY[i],
                       col=colour[i])
        if (i<ns) graphics::lines(rep(X[i+1],2),
                                  c(Y[i]-fact*sY[i],Y[i+1]+fact*sY[i+1]))
    }
    if (plateau){
        plat$n <- nrow(x)
        if (title) graphics::title(plateau.title(plat,sigdig=sigdig,Ar=FALSE))
        return(invisible(plat))
    }
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
                             plateau.col=rgb(0,1,0,0.5),
                             non.plateau.col=rgb(0,1,1,0.5),sigdig=2,
                             exterr=TRUE,line.col='red',lwd=2,
                             i2i=FALSE,...){
    tt <- ArAr.age(x,jcu=FALSE,exterr=FALSE,i2i=i2i)
    X <- cbind(x$x[,'Ar39'],tt)
    x.lab <- expression(paste("cumulative ",""^"39","Ar fraction"))
    plat <- agespectrum.default(X,alpha=alpha,xlab=x.lab,ylab='age [Ma]',
                                plateau=plateau,sigdig=sigdig,
                                line.col=line.col,
                                lwd=lwd,title=FALSE,...)
    if (plateau){
        out <- plat
        # calculate the weighted mean Ar40Ar39 ratio from the weighted mean age
        R <- get.ArAr.ratio(plat$mean['x'],plat$mean['s[x]'],x$J[1],0,exterr=FALSE)
        # recalculate the weighted mean age, this time
        # taking into account decay and J uncertainties
        out$mean[1:2] <- get.ArAr.age(R[1],R[2],x$J[1],x$J[2],exterr=exterr)
        out$mean[3] <- nfact(alpha)*out$mean[2]
        graphics::title(plateau.title(out,sigdig=sigdig,Ar=TRUE,units='Ma'))
        return(invisible(out))
    }
}

# x is a three column vector with Ar39
# cumulative fractions, ages and uncertainties
plateau <- function(x,alpha=0.05){
    X <- x[,1]/sum(x[,1],na.rm=TRUE)
    YsY <- subset(x,select=c(2,3))
    ns <- length(X)
    out <- list()
    out$mean <- YsY[1,1:2]
    out$mswd <- 1
    out$p.value <- 1
    out$fract <- 0
    for (i in 1:(ns-1)){ # at least two steps in a plateau
        for (j in (i+1):ns){
            fract <- sum(X[i:j],na.rm=TRUE)
            Y <- YsY[i:j,1]
            sY <- YsY[i:j,2]
            valid <- chauvenet(Y,sY,valid=rep(TRUE,j-i+1))
            if (any(!valid)){
                break;
            } else if (fract > out$fract){
                out <- weightedmean(YsY[i:j,],plot=FALSE,
                                    detect.outliers=FALSE)
                out$i <- i:j
                out$fract <- fract
            }
        }
    }
    out
}

plateau.title <- function(fit,sigdig=2,Ar=TRUE,units=''){
    rounded.mean <- roundit(fit$mean[1],fit$mean[2:3],sigdig=sigdig)
    line1 <- substitute('mean ='~a%+-%b~'|'~c~u~'(n='~n/N~')',
                        list(a=rounded.mean[1],
                             b=rounded.mean[2],
                             c=rounded.mean[3],
                             u=units,
                             n=length(fit$i),
                             N=fit$n))
    a <- signif(100*fit$fract,sigdig)
    if (Ar)
        line2 <- bquote(paste("Includes ",.(a),"% of the ",""^"39","Ar"))
    else
        line2 <- bquote(paste("Includes ",.(a),"% of the spectrum"))
    graphics::mtext(line1,line=1)
    graphics::mtext(line2,line=0)
}
