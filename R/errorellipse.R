#' Get coordinates of error ellipse for plotting
#'
#' Construct an error ellipse age a given confidence level from its
#' centre and covariance matrix
#'
#' @param x x-coordinate (scalar) for the centre of the ellipse
#' @param y y-coordinate (scalar) for the centre of the ellipse
#' @param covmat covariance matrix of the x-y coordinates
#' @param alpha the probability cutoff for the error ellipses
#' @param n the resolution of the error ellipses
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

# x = matrix with columns X, sX, Y, sY, rXY
scatterplot <- function(x,xlim=NA,ylim=NA,alpha=0.05,
                        show.numbers=FALSE,show.ellipses=1,levels=NA,
                        ellipse.col=c("#00FF0080","#FF000080"),a=NA,
                        b=NA,line.col='red',lwd=2, new.plot=TRUE,
                        empty=FALSE,...){
    colnames(x) <- c('X','sX','Y','sY','rXY')
    if (any(is.na(xlim))) xlim <- get.limits(x[,'X'],x[,'sX'])
    if (any(is.na(ylim))) ylim <- get.limits(x[,'Y'],x[,'sY'])
    if (new.plot) graphics::plot(xlim,ylim,type='n',xlab='',ylab='')
    if (new.plot & empty) return()
    if (!is.na(a) & !is.na(b)){
        graphics::lines(xlim,a+b*xlim,col=line.col,lwd=lwd)
    }
    ns <- nrow(x)
    x0 <- x[,'X']
    y0 <- x[,'Y']
    if (show.ellipses==0){
        if (show.numbers) graphics::text(x0,y0,1:ns,...)
        else graphics::points(x0,y0,...)
    } else if (show.ellipses==1){
        ellipse.cols <- set.ellipse.colours(ns=ns,levels=levels,col=ellipse.col)
        for (i in 1:ns){
            if (!any(is.na(x[i,]))){
                covmat <- cor2cov2(x[i,'sX'],x[i,'sY'],x[i,'rXY'])
                ell <- ellipse(x0[i],y0[i],covmat,alpha=alpha)
                graphics::polygon(ell,col=ellipse.cols[i])
                if (show.numbers) graphics::text(x0[i],y0[i],i)
                else graphics::points(x0[i],y0[i],pch=19,cex=0.25)
            }
        }
        colourbar(z=levels,col=ellipse.col)
    } else {
        if (show.numbers) graphics::text(x0,y0,1:ns,adj=c(0,1))
        else graphics::points(x0,y0,pch=19,cex=0.5)
        fact <- stats::qnorm(1-alpha/2)
        dx <- fact*x[,'sX']
        dy <- fact*x[,'sY']
        graphics::arrows(x0,y0-dy,x0,y0+dy,code=3,angle=90,length=0.05)
        graphics::arrows(x0-dx,y0,x0+dx,y0,code=3,angle=90,length=0.05)
    }
}
