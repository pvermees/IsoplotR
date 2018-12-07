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

# d = matrix with columns X, sX, Y, sY, rXY
scatterplot <- function(d,xlim=NA,ylim=NA,alpha=0.05,
                        show.numbers=FALSE,show.ellipses=1,levels=NA,
                        clabel="",ellipse.col=c("#00FF0080","#FF000080"),
                        fit='none',new.plot=TRUE,ci.col='gray80',
                        line.col='black',lwd=1,empty=FALSE,
                        hide=NULL,omit=NULL,omit.col=NA,
                        addcolourbar=TRUE,...){
    ns <- nrow(d)
    sn <- 1:ns
    plotit <- sn%ni%hide
    calcit <- sn%ni%c(hide,omit)
    colnames(d) <- c('X','sX','Y','sY','rXY')
    if (any(is.na(xlim))) xlim <- get.limits(d[plotit,'X'],d[plotit,'sX'])
    if (any(is.na(ylim))) ylim <- get.limits(d[plotit,'Y'],d[plotit,'sY'])
    if (new.plot) graphics::plot(xlim,ylim,type='n',xlab='',ylab='',bty='n')
    if (new.plot & empty) return()
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
        plot_points(d[,'X'],d[,'Y'],mybg=colour,mycex=1,
                    show.numbers=show.numbers,hide=hide,omit=omit,...)
    } else if (show.ellipses==1){ # error ellipse
        for (i in sn[plotit]){
            if (!any(is.na(d[i,]))){
                covmat <- cor2cov2(d[i,'sX'],d[i,'sY'],d[i,'rXY'])
                ell <- ellipse(d[i,'X'],d[i,'Y'],covmat,alpha=alpha)
                graphics::polygon(ell,col=colour[i])
                if (show.numbers) graphics::text(d[i,'X'],d[i,'Y'],i)
                else graphics::points(d[i,'X'],d[i,'Y'],pch=19,cex=0.25)
            }
        }
    } else { # error cross
        if (show.numbers)
            graphics::text(d[plotit,'X'],d[plotit,'Y'],
                           sn[plotit],adj=c(0,1),col=colour)
        else
            graphics::points(d[plotit,'X'],d[plotit,'Y'],
                             pch=19,cex=0.5,col=colour)
        fact <- stats::qnorm(1-alpha/2)
        dx <- fact*d[plotit,'sX']
        dy <- fact*d[plotit,'sY']
        graphics::arrows(d[plotit,'X'],d[plotit,'Y']-dy,
                         d[plotit,'X'],d[plotit,'Y']+dy,
                         code=3,angle=90,
                         length=0.05,col=colour)
        graphics::arrows(d[plotit,'X']-dx,d[plotit,'Y'],
                         d[plotit,'X']+dx,d[plotit,'Y'],
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
