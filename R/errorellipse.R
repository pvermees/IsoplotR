#' @title Get error ellipse coordinates for plotting
#'
#' @description Constructs an error ellipse at a given confidence level
#' from its centre and covariance matrix
#'
#' @param x x-coordinate (scalar) for the centre of the ellipse
#' @param y y-coordinate (scalar) for the centre of the ellipse
#' @param covmat the [\code{2x2}] covariance matrix of the x-y coordinates
#' @param alpha the probability cutoff for the error ellipses
#' @param n the resolution (number of segments) of the error ellipses
#'
#' @return an [\code{nx2}] matrix of plot coordinates
#' @examples
#' x = 99; y = 101;
#' covmat <- matrix(c(1,0.9,0.9,1),nrow=2)
#' ell <- ellipse(x,y,covmat)
#' plot(c(90,110),c(90,110),type='l')
#' polygon(ell,col=rgb(0,1,0,0.5))
#' points(x,y,pch=21,bg='black')
#' @export
ellipse <- function(x,y,covmat,alpha=0.05,n=50){
    cutoff <- stats::qchisq(1-alpha,2)
    e <- eigen(covmat)
    a <- sqrt(cutoff*abs(e$values[1])) # major axis
    b <- sqrt(cutoff*abs(e$values[2])) # minor axis
    v <- e$vectors[,1] # largest eigenvector
    beta <- atan(v[2]/v[1]) # rotation angle of the ellipse
    theta <- seq(0, 2 * pi, length=n)
    out <- matrix(0,nrow=n,ncol=2)
    out[,1] <- x + a * cos(theta) * cos(beta) - b * sin(theta) * sin(beta)
    out[,2] <- y + a * cos(theta) * sin(beta) + b * sin(theta) * cos(beta)
    colnames(out) <- c('x','y')
    out
}

#' @title Create a scatter plot with error ellipses or crosses
#'
#' @description
#' Takes bivariate data with (correlated) uncertainties as input and
#' produces a scatter plot with error ellipses or crosses as output.
#' (optionally) displays the linear fit on this diagram, and can show
#' a third variable as a colour scale.
#'
#' @param xy matrix with columns \code{X, sX, Y, sY(, rXY)}
#' @param alpha the probability cutoff for the error ellipses
#' @param show.numbers logical flag (\code{TRUE} to show grain
#'     numbers)
#' @param show.ellipses show the data as:
#' 
#' \code{0}: points
#' 
#' \code{1}: error ellipses
#' 
#' \code{2}: error crosses
#' 
#' @param levels a vector with additional values to be displayed as
#'     different background colours within the error ellipses.
#' @param clabel label for the colour scale
#' @param ellipse.fill
#' Fill colour for the error ellipses. This can either be a single
#' colour or multiple colours to form a colour ramp. Examples:
#'
#' a single colour: \code{rgb(0,1,0,0.5)}, \code{'#FF000080'},
#' \code{'white'}, etc.;
#'
#' multiple colours: \code{c(rbg(1,0,0,0.5)},
#' \code{rgb(0,1,0,0.5))}, \code{c('#FF000080','#00FF0080')},
#' \code{c('blue','red')}, \code{c('blue','yellow','red')}, etc.;
#'
#' a colour palette: \code{rainbow(n=100)},
#' \code{topo.colors(n=100,alpha=0.5)}, etc.; or
#'
#' a reversed palette: \code{rev(topo.colors(n=100,alpha=0.5))},
#' etc.
#'
#' For empty ellipses, set \code{ellipse.col=NA}
#' @param ellipse.stroke the stroke colour for the error
#'     ellipses. Follows the same formatting guidelines as
#'     \code{ellipse.fill}
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
#' @param omit.fill fill colour that should be used for the omitted
#'     aliquots.
#' @param omit.stroke stroke colour that should be used for the
#'     omitted aliquots.
#' @param addcolourbar add a colour bar to display the colours used to
#'     \code{levels}
#' @param bg background colour for the plot symbols (only used if
#'     \code{show.ellipses=0}).
#' @param cex plot symbol magnification.
#' @param xlim (optional) two-element vector with the x-axis limits
#' @param ylim (optional) two-element vector with the y-axis limits
#' @param xlab (optional) x-axis label (only used when
#'     \code{add=FALSE})
#' @param ylab (optional) y-axis label (only used when
#'     \code{add=FALSE})
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
scatterplot <- function(xy,alpha=0.05,show.numbers=FALSE,
                        show.ellipses=1,levels=NA,clabel="",
                        ellipse.fill=c("#00FF0080","#FF000080"),
                        ellipse.stroke="black",fit='none',add=FALSE,
                        empty=FALSE, ci.col='gray80',line.col='black',
                        lwd=1,hide=NULL,omit=NULL,omit.fill=NA,
                        omit.stroke="grey",addcolourbar=TRUE,
                        bg,cex,xlim=NULL,ylim=NULL,xlab,ylab,...){
    ns <- nrow(xy)
    if (ncol(xy)==4) xy <- cbind(xy,rep(0,ns))
    sn <- 1:ns
    plotit <- sn%ni%hide
    calcit <- sn%ni%c(hide,omit)
    colnames(xy) <- c('X','sX','Y','sY','rXY')
    if (is.null(xlim)) xlim <- get.limits(xy[plotit,'X'],xy[plotit,'sX'])
    if (is.null(ylim)) ylim <- get.limits(xy[plotit,'Y'],xy[plotit,'sY'])
    if (!add){
        if (missing(xlab)) xlab <- ''
        if (missing(ylab)) ylab <- ''
        graphics::plot(xlim,ylim,type='n',xlab=xlab,ylab=ylab,bty='n')
        if (empty) return()
    }
    if (!identical(fit,'none')){
        plot_isochron_line(fit,x=seq(xlim[1],xlim[2],length.out=100),
                           ci.col=ci.col,col=line.col,lwd=lwd)
    }
    graphics::box()
    nolevels <- all(is.na(levels))
    fill <- set.ellipse.colours(ns=ns,levels=levels,col=ellipse.fill,
                                hide=hide,omit=omit,omit.col=omit.fill)
    stroke <- set.ellipse.colours(ns=ns,levels=levels,col=ellipse.stroke,
                                  hide=hide,omit=omit,omit.col=omit.stroke)
    if (show.ellipses==0){ # points and or text
        if (missing(cex)) cex <- 1
        if (missing(bg)) bg <- fill
        plot_points(xy[,'X'],xy[,'Y'],bg=bg,cex=cex,col=stroke,
                    show.numbers=show.numbers,hide=hide,omit=omit,...)
    } else if (show.ellipses==1){ # error ellipse
        if (missing(cex)) cex <- 0.25
        for (i in sn[plotit]){
            if (!any(is.na(xy[i,]))){
                covmat <- cor2cov2(xy[i,'sX'],xy[i,'sY'],xy[i,'rXY'])
                ell <- ellipse(xy[i,'X'],xy[i,'Y'],covmat,alpha=alpha)
                graphics::polygon(ell,col=fill[i],border=stroke[i])
                if (show.numbers) graphics::text(xy[i,'X'],xy[i,'Y'],i)
                else graphics::points(xy[i,'X'],xy[i,'Y'],pch=19,cex=cex)
            }
        }
    } else { # error cross
        if (missing(cex)) cex <- 0.5
        if (show.numbers)
            graphics::text(xy[plotit,'X'],xy[plotit,'Y'],
                           sn[plotit],adj=c(0,1),col=stroke)
        else
            graphics::points(xy[plotit,'X'],xy[plotit,'Y'],
                             pch=21,cex=cex,col=stroke,bg=fill)
        fact <- stats::qnorm(1-alpha/2)
        dx <- fact*xy[plotit,'sX']
        dy <- fact*xy[plotit,'sY']
        graphics::arrows(xy[plotit,'X'],xy[plotit,'Y']-dy,
                         xy[plotit,'X'],xy[plotit,'Y']+dy,
                         code=3,angle=90,length=0.05,col=stroke)
        graphics::arrows(xy[plotit,'X']-dx,xy[plotit,'Y'],
                         xy[plotit,'X']+dx,xy[plotit,'Y'],
                         code=3,angle=90,length=0.05,col=stroke)
    }
    if (!nolevels & addcolourbar){
        colourbar(z=levels[calcit],fill=ellipse.fill,
                  stroke=ellipse.stroke,clabel=clabel)
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
