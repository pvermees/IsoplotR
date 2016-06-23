isochron <- function(x,...){ UseMethod("isochron",x) }
# x is a list with the following vectors:
# X the x-variable
# Y the y-variable
# sX the standard error of X
# sY the standard error of Y
# rXY the correlation coefficient of X and Y
isochron.default <- function(x,xlim=NA,ylim=NA,alpha=0.05,
                             show.numbers=FALSE,ellipse.col=rgb(0,1,0,0.5),...){
    fit <- yorkfit(x$X,x$Y,x$sX,x$sY,x$rXY)
    if (is.na(xlim)) xlim <- get.limits(x$X,x$sX)
    if (is.na(ylim)) ylim <- get.limits(x$Y,x$sY)
    plot(xlim,ylim,type='n')
    lines(xlim,fit$a+fit$b*xlim)
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
    title(paste0('MSWD=',fit$mswd))
}
isochron.ArAr <- function(x,xlim=NA,ylim=NA,alpha=0.05,
                          show.numbers=FALSE,ellipse.col=rgb(0,1,0,0.5),
                          inverse=TRUE,...){
    d <- data2york(x,get.selection(x,inverse))
    fit <- yorkfit(d$X,d$Y,d$sX,d$sY,d$rXY)
    isochron.default(d,xlim,ylim,alpha,show.numbers,
                     ellipse.col,wetherill,...)
}

get.limits <- function(X,sX){
    minx <- min(X-3*sX)
    maxx <- max(X+3*sX)    
    c(minx,maxx)
}

isochron.title <- function(fit){
    rounded.age <- roundit(fit$age,fit$age.err)
    line1 <- substitute('concordia age ='~a%+-%b~'[Ma] (1'~sigma~')',
                        list(a=rounded.age$x, b=rounded.age$err))
    line2 <- substitute('MSWD (concordance) ='~a~', p('~chi^2*')='~b,
                        list(a=signif(fit$mswd$concordance,2),
                             b=signif(fit$p.value$concordance,2)))
    graphics::mtext(line1,line=1)
    graphics::mtext(line2,line=0)
}
