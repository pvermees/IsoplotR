#' Calculate and plot isochrons
#'
#' Plots cogenetic \eqn{^{40}}Ar/\eqn{^{39}}Ar data as X-Y
#' scatterplots, fits an isochron curve through them using the
#' \code{yorkfit} function, and computes the corresponding isochron
#' age, including decay constant uncertainties.
#'
#' @param x EITHER a list or a matrix with the following vectors:
#'
#' \code{X:} the x-variable
#'
#' \code{Y:} the y-variable
#'
#' \code{sX:} the standard error of X
#'
#' \code{sY:} the standard error of Y
#'
#' \code{rXY:} the correlation coefficient of X and Y
#'
#' OR an object of class \code{ArAr}
#'
#' @param xlim 2-element vector with the plot limits of the x-axis
#' @param ylim 2-element vector with the plot limits of the y-axis
#' @param alpha confidence cutoff for the error ellipses
#' @param show.numbers boolean flag (TRUE to show grain numbers)
#' @param sigdig the number of significant digits of the numerical
#'     values reported in the title of the graphical output
#' @param ellipse.col background colour of the error ellipses
#' @param line.col colour of the isochron line
#' @param lwd line width
#' @param title add a title to the plot?
#' @param ... optional arguments
#' @rdname isochron
#' @export
isochron <- function(x,...){ UseMethod("isochron",x) }
#' @rdname isochron
#' @export
isochron.default <- function(x,xlim=NA,ylim=NA,alpha=0.05,
                             sigdig=2,show.numbers=FALSE,
                             ellipse.col=rgb(0,1,0,0.5),
                             line.col='red',lwd=2,title=TRUE,...){
    if (hasClass(x,'matrix') | hasClass(x,'data.frame'))
        x <- list(X=x[,1],sX=x[,2],Y=x[,3],sY=x[,4],rXY=x[,5])
    fit <- yorkfit(x$X,x$Y,x$sX,x$sY,x$rXY)
    scatterplot(x,xlim=xlim,ylim=ylim,alpha=alpha,
                show.numbers=show.numbers, ellipse.col=ellipse.col,
                a=fit$a[1],b=fit$b[1], line.col=line.col,lwd=lwd)
    if (title)
        title(regression.title(fit,sigdig=sigdig),xlab='X',ylab='Y')
}
#' @param plot if \code{FALSE}, suppresses the graphical output
#' @param inverse if \code{TRUE}, plots \eqn{^{36}}Ar/\eqn{^{40}}Ar
#'     vs. \eqn{^{39}}Ar/\eqn{^{40}}Ar. If \code{FALSE}, plots
#'     \eqn{^{40}}Ar/\eqn{^{36}}Ar vs. \eqn{^{39}}Ar/\eqn{^{36}}Ar.
#' @return
#' if \code{plot=FALSE}, returns a list with the following items:
#'
#' \code{a:} the intercept of the straight line fit and its standard error
#' 
#' \code{b:} the slope of the fit and its standard error
#' 
#' \code{y0:} the atmospheric \eqn{^{40}}Ar/\eqn{^{36}}Ar ratio and its standard error
#' 
#' \code{age:} the \eqn{^{40}}Ar/\eqn{^{39}}Ar age and its standard error
#' 
#' @examples
#' data(examples)
#' isochron(examples$ArAr)
#' @rdname isochron
#' @export
isochron.ArAr <- function(x,xlim=NA,ylim=NA,alpha=0.05,sigdig=sigdig,
                          show.numbers=FALSE,ellipse.col=rgb(0,1,0,0.5),
                          inverse=TRUE,line.col='red',lwd=2,plot=TRUE,...){
    d <- data2york(x,get.selection(x,inverse))
    if (inverse){
        fit <- yorkfit(d$Y,d$X,d$sY,d$sX,-d$rXY) # X and Y reversed!
        x0 <- 1/fit$a[1]
        sx0 <- fit$a[2]/fit$a[1]^2
        y0 <- -fit$b[1]/fit$a[1]
        sy0 <- y0*sqrt((fit$a[2]/fit$a[1])^2 + (fit$b[2]/fit$b[1])^2)
        tt <- get.ArAr.age(x0,sx0,x$J[1],x$J[2])
        x.lab <- expression(paste(""^"39","Ar/"^"40","Ar"))
        y.lab <- expression(paste(""^"36","Ar/"^"40","Ar"))
    } else {
        fit <- yorkfit(d$X,d$Y,d$sX,d$sY,d$rXY)
        y0 <- fit$a[1]
        sy0 <- fit$a[2]
        tt <- get.ArAr.age(fit$b[1],fit$b[2],x$J[1],x$J[2])
        x.lab <- expression(paste(""^"39","Ar/"^"36","Ar"))
        y.lab <- expression(paste(""^"40","Ar/"^"36","Ar"))
    }
    out <- fit
    class(out) <- "isochron"
    out$y0 <- c(y0,sy0)
    out$age <- tt
    if (plot) {
        isochron.default(d,alpha=alpha,show.numbers=show.numbers,
                         ellipse.col=ellipse.col,a=fit$a[1],b=fit$b[1],
                         line.col=line.col,lwd=lwd,title=FALSE)
        tt <- roundit(out$age[1],out$age[2])
        title(isochron.title(out,sigdig=2),xlab=x.lab,ylab=y.lab)
    } else {
        return(out)
    }
}

get.limits <- function(X,sX){
    minx <- min(X-3*sX)
    maxx <- max(X+3*sX)    
    c(minx,maxx)
}

isochron.title <- function(fit,sigdig=2){
    rounded.age <- roundit(fit$age[1],fit$age[2],sigdig=sigdig)
    rounded.intercept <- roundit(fit$y0[1],fit$y0[2],sigdig=sigdig)
    line1 <- substitute('age ='~a%+-%b~'(1'~sigma~'), intercept ='~c%+-%d~'(1'~sigma~')',
                        list(a=rounded.age$x, b=rounded.age$err,
                             c=rounded.intercept$x, d=rounded.intercept$err))
    line2 <- substitute('MSWD ='~a~', p('~chi^2*')='~b,
                        list(a=signif(fit$mswd,2), b=signif(fit$p.value,2)))
    graphics::mtext(line1,line=1)
    graphics::mtext(line2,line=0)
}

regression.title <- function(fit,sigdig=2){
    intercept <- roundit(fit$a[1],fit$a[2],sigdig=sigdig)
    slope <- roundit(fit$b[1],fit$b[2],sigdig=sigdig)
    line1 <- substitute('slope ='~a%+-%b~'(1'~sigma~'), intercept ='~c%+-%d~'(1'~sigma~')',
                        list(a=slope$x, b=slope$err,
                             c=intercept$x, d=intercept$err))
    line2 <- substitute('MSWD ='~a~', p('~chi^2*')='~b,
                        list(a=signif(fit$mswd,2), b=signif(fit$p.value,2)))
    graphics::mtext(line1,line=1)
    graphics::mtext(line2,line=0)
}
