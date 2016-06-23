isochron <- function(x,...){ UseMethod("isochron",x) }
# x is a list with the following vectors:
# X the x-variable
# Y the y-variable
# sX the standard error of X
# sY the standard error of Y
# rXY the correlation coefficient of X and Y
isochron.default <- function(x,plot=TRUE,xlim=NA,ylim=NA, alpha=0.05,
                             show.numbers=FALSE,
                             ellipse.col=rgb(0,1,0,0.5),...){
    fit <- yorkfit(x$X,x$Y,x$sX,x$sY,x$rXY)
    if (is.na(xlim)) xlim <- get.limits(x$X,x$sX)
    if (is.na(ylim)) ylim <- get.limits(x$Y,x$sY)
    plot(xlim,ylim,type='n')
    ns <- length(x$X)
    for (i in 1:ns){
        x0 <- x$X[i]
        y0 <- x$Y[i]
        covmat <- cor2cov(x$sX[i],x$sY[i],x$rXY[i])
        ell <- ellipse(x0,y0,covmat,alpha=alpha)
        polygon(ell,col=ellipse.col)
        points(x0,y0,pch=19,cex=0.25)
        if (show.numbers) { text(x0,y0,i) }
    }
}
isochron.UPb <- function(x,plot=TRUE,wetherill=TRUE,...){
    d <- UPb2york(x,wetherill=wetherill)
    isochron.default(d,plot=plot,...)
}

get.limits <- function(X,sX){
    minx <- min(X-3*sX)
    maxx <- max(X+3*sX)
    c(minx,maxx)
}
