isochron <- function(x,...){ UseMethod("isochron",x) }
# x is a list with the following vectors:
# X the x-variable
# Y the y-variable
# sX the standard error of X
# sY the standard error of Y
# rXY the correlation coefficient of X and Y
isochron.default <- function(x,xlim=NA,ylim=NA,alpha=0.05,
                             show.numbers=FALSE,
                             ellipse.col=rgb(0,1,0,0.5),
                             line.col='grey',lwd=2,...){
    fit <- yorkfit(x$X,x$Y,x$sX,x$sY,x$rXY)
    if (any(is.na(xlim))) xlim <- get.limits(x$X,x$sX)
    if (any(is.na(ylim))) ylim <- get.limits(x$Y,x$sY)
    plot(xlim,ylim,type='n',xlab='',ylab='')
    lines(xlim,fit$a[1]+fit$b[1]*xlim,col=line.col,lwd=lwd)
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
isochron.ArAr <- function(x,xlim=NA,ylim=NA,alpha=0.05,
                          show.numbers=FALSE,ellipse.col=rgb(0,1,0,0.5),
                          inverse=TRUE,plot=TRUE,...){
    d <- data2york(x,get.selection(x,inverse))
    fit <- yorkfit(d$X,d$Y,d$sX,d$sY,d$rXY)
    out <- fit
    class(out) <- "isochron"
    out$inverse <- inverse
    if (inverse){
        x0 <- -fit$a[1]/fit$b[1]
        # doesn't take into account covariance between a and b!
        sx0 <- x0*sqrt((fit$a[2]/fit$a[1])^2 + (fit$b[2]/fit$b[1])^2)
        y0 <- 1/fit$a[1]
        sy0 <- fit$a[2]/fit$a[1]^2
        out$y0 <- c(y0,sy0)
        out$age <- get.ArAr.age(1/x0,sx0/x0^2,x$J[1],x$J[2])
        x.lab <- expression(paste(""^"39","Ar/"^"40","Ar"))
        y.lab <- expression(paste(""^"36","Ar/"^"40","Ar"))
    } else {
        out$y0 <- fit$a
        out$age <- get.ArAr.age(fit$b[1],fit$b[2],x$J[1],x$J[2])
        x.lab <- expression(paste(""^"39","Ar/"^"36","Ar"))
        y.lab <- expression(paste(""^"40","Ar/"^"36","Ar"))
    }
    if (plot) {
        isochron.default(d,xlim,ylim,alpha,show.numbers,
                         ellipse.col,...)
        tt <- roundit(out$age[1],out$age[2])
        title(isochron.title(out),xlab=x.lab,ylab=y.lab)
    } else {
        return(out)
    }
}

get.limits <- function(X,sX){
    minx <- min(X-3*sX)
    maxx <- max(X+3*sX)    
    c(minx,maxx)
}

isochron.title <- function(fit){
    rounded.age <- roundit(fit$age[1],fit$age[2])
    rounded.intercept <- roundit(fit$y0[1],fit$y0[2])
    line1 <- substitute('age ='~a%+-%b~'(1'~sigma~'), intercept ='~c%+-%d~'(1'~sigma~')',
                        list(a=rounded.age$x, b=rounded.age$err,
                             c=rounded.intercept$x, d=rounded.intercept$err))
    line2 <- substitute('MSWD ='~a~', p('~chi^2*')='~b,
                        list(a=signif(fit$mswd,2), b=signif(fit$p.value,2)))
    graphics::mtext(line1,line=1)
    graphics::mtext(line2,line=0)
}
