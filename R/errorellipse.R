#' Get coordinates of error ellipse for plotting
#'
#' Constructs an error ellipse at a given confidence level from its
#' centre and covariance matrix
#'
#' @param x x-coordinate (scalar) for the centre of the ellipse
#' @param y y-coordinate (scalar) for the centre of the ellipse
#' @param covmat the [\code{2 x 2}] covariance matrix of the x-y coordinates
#' @param alpha the probability cutoff for the error ellipses
#' @param n the resolution (number of segments) of the error ellipses
#'
#' @return an [\code{n x 2}] matrix of plot coordinates
#' @examples
#' x = 99; y = 101;
#' covmat <- matrix(c(1,0.9,0.9,1),nrow=2)
#' ell <- ellipse(x,y,covmat)
#' plot(c(90,110),c(90,110),type='l')
#' polygon(ell,col=rgb(0,1,0,0.5))
#' points(x,y,pch=21,bg='black')
#' @export
ellipse <- function(x,y,covmat,alpha=0.05,n=50){
    nn <- 50
    cutoff <- stats::qchisq(1-alpha,2)
    e <- eigen(covmat)
    a <- sqrt(cutoff*abs(e$values[1])) # major axis
    b <- sqrt(cutoff*abs(e$values[2])) # minor axis
    v <- e$vectors[,1] # largest eigenvector
    beta <- atan(v[2]/v[1]) # rotation angle of the ellipse
    theta <- seq(0, 2 * pi, length=nn)
    out <- matrix(0,nrow=nn,ncol=2)
    out[,1] <- x + a * cos(theta) * cos(beta) - b * sin(theta) * sin(beta)
    out[,2] <- y + a * cos(theta) * sin(beta) + b * sin(theta) * cos(beta)
    colnames(out) <- c('x','y')
    out
}

#' Create a scatter plot with error ellipses or crosses
#'
#' Takes bivariate data with (correlated) uncertainties as input and
#' produces a scatter plot with error ellipses or crosses as output.
#' (optionally) displays the linear fit on this diagram, and can show
#' a third variable as a colour scale.
#'
#' @param xy matrix with columns X, sX, Y, sY(, rXY)
#' @param xlim (optional) two-element vector with the x-axis limits
#' @param ylim (optional) two-element vector with the y-axis limits
#' @param alpha the probability cutoff for the error ellipses
#' @param show.numbers logical flag (\code{TRUE} to show grain
#'     numbers)
#' @param show.ellipses show the data as: \enumerate{ \item points
#'     \item error ellipses \item error crosses }
#' @param levels a vector with additional values to be displayed as
#'     different background colours within the error ellipses.
#' @param clabel label for the colour scale
#' @param ellipse.col a vector of two background colours for the error
#'     ellipses. If \code{levels=NA}, then only the first colour will
#'     be used. If \code{levels} is a vector of numbers, then
#'     \code{ellipse.col} is used to construct a colour ramp.
#' @param fit the output of \code{york()} (optional).
#' @param add if \code{TRUE}, adds the points and lines to the
#'     existing plot.
#' @param empty set up an empty plot with the right axis limits to fit
#'     the data
#' @param ci.col the fill colour for the confidence interval of the
#'     intercept and slope.
#' @param line.col colour of the regression line
#' @param lwd line width of the regression line
#' @param hide vector with indices of aliquots that should be removed
#'     from the plot.
#' @param omit vector with indices of aliquots that should be plotted
#'     but omitted from the isochron age calculation.
#' @param omit.col colour that should be used for the omitted
#'     aliquots.
#' @param addcolourbar add a colour bar to display the colours used to
#'     \code{levels}
#' @param ... optional arguments to format the points and text.
#' 
#' @examples
#' X <- c(1.550,12.395,20.445,20.435,20.610,24.900,
#'        28.530,50.540,51.595,86.51,106.40,157.35)
#' Y <- c(.7268,.7809,.8200,.8116,.8160,.8302,
#'        .8642,.9534,.9617,1.105,1.230,1.440)
#' sX <- X*0.02
#' sY <- Y*0.01
#' dat <- cbind(X,sX,Y,sY)
#' scatterplot(dat,fit=york(dat),show.ellipses=2)
#' @export
scatterplot <- function(xy,xlim=NA,ylim=NA,alpha=0.05,
                        show.numbers=FALSE,show.ellipses=1,levels=NA,
                        clabel="",ellipse.col=c("#00FF0080","#FF000080"),
                        fit='none',add=FALSE,empty=FALSE,
                        ci.col='gray80',line.col='black',lwd=1,
                        hide=NULL,omit=NULL,omit.col=NA,
                        addcolourbar=TRUE,...){
    ns <- nrow(xy)
    if (ncol(xy)==4) xy <- cbind(xy,rep(0,ns))
    sn <- 1:ns
    plotit <- sn%ni%hide
    calcit <- sn%ni%c(hide,omit)
    colnames(xy) <- c('X','sX','Y','sY','rXY')
    if (any(is.na(xlim))) xlim <- get.limits(xy[plotit,'X'],xy[plotit,'sX'])
    if (any(is.na(ylim))) ylim <- get.limits(xy[plotit,'Y'],xy[plotit,'sY'])
    if (!add) graphics::plot(xlim,ylim,type='n',xlab='',ylab='',bty='n')
    if (!add & empty) return()
    if (!identical(fit,'none'))
        plot_isochron_line(fit,x=seq(xlim[1],xlim[2],length.out=100),
                           ci.col=ci.col,col=line.col,lwd=lwd)
    graphics::box()
    haslevels <- !all(is.na(levels))
    if (show.ellipses==2 && all(is.na(levels))){
        colour <- rep('black',ns)
    } else {
        colour <- set.ellipse.colours(ns=ns,levels=levels,col=ellipse.col,
                                      hide=hide,omit=omit,omit.col=omit.col)
    }
    if (show.ellipses==0){ # points and or text
        plot_points(xy[,'X'],xy[,'Y'],mybg=colour,mycex=1,
                    show.numbers=show.numbers,hide=hide,omit=omit,...)
    } else if (show.ellipses==1){ # error ellipse
        for (i in sn[plotit]){
            if (!any(is.na(xy[i,]))){
                covmat <- cor2cov2(xy[i,'sX'],xy[i,'sY'],xy[i,'rXY'])
                ell <- ellipse(xy[i,'X'],xy[i,'Y'],covmat,alpha=alpha)
                graphics::polygon(ell,col=colour[i])
                if (show.numbers) graphics::text(xy[i,'X'],xy[i,'Y'],i)
                else graphics::points(xy[i,'X'],xy[i,'Y'],pch=19,cex=0.25)
            }
        }
    } else { # error cross
        if (show.numbers)
            graphics::text(xy[plotit,'X'],xy[plotit,'Y'],
                           sn[plotit],adj=c(0,1),col=colour)
        else
            graphics::points(xy[plotit,'X'],xy[plotit,'Y'],
                             pch=19,cex=0.5,col=colour)
        fact <- stats::qnorm(1-alpha/2)
        dx <- fact*xy[plotit,'sX']
        dy <- fact*xy[plotit,'sY']
        graphics::arrows(xy[plotit,'X'],xy[plotit,'Y']-dy,
                         xy[plotit,'X'],xy[plotit,'Y']+dy,
                         code=3,angle=90,
                         length=0.05,col=colour)
        graphics::arrows(xy[plotit,'X']-dx,xy[plotit,'Y'],
                         xy[plotit,'X']+dx,xy[plotit,'Y'],
                         code=3,angle=90,
                         length=0.05,col=colour)
    }
    if (haslevels & addcolourbar){
        colourbar(z=levels[calcit],col=ellipse.col,clabel=clabel)
    }
}

plot_isochron_line <- function(fit,x,ci.col='gray80',...){
    y <- fit$a[1]+fit$b[1]*x
    e <- fit$fact*sqrt(fit$a[2]^2 + 2*x*fit$cov.ab + (fit$b[2]*x)^2)
    cix <- c(x,rev(x))
    ciy <- c(y+e,rev(y-e))
    graphics::polygon(cix,ciy,col=ci.col,border=NA)
    graphics::lines(x,y,...)
}
