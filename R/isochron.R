#' Calculate and plot isochrons
#'
#' Plots cogenetic \eqn{^{40}}Ar/\eqn{^{39}}Ar data as X-Y
#' scatterplots, fits an isochron curve through them using the
#' \code{yorkfit} function, and computes the corresponding isochron
#' age, including decay constant uncertainties.
#'
#' @param x EITHER a list or a matrix with the following vectors or
#'     columns:
#' \describe{
#' \item{X}{the x-variable}
#' \item{sX}{the standard error of \code{X}}
#' \item{Y}{the y-variable}
#' \item{sY}{the standard error of \code{Y}}
#' \item{rXY}{the correlation coefficient of \code{X} and \code{Y}}
#' }
#' OR
#'
#' an object of class \code{ArAr}, \code{ReOs}, \code{RbSr} or
#' \code{SmNd}
#'
#' @param xlim 2-element vector with the plot limits of the x-axis
#' @param ylim 2-element vector with the plot limits of the y-axis
#' @param alpha confidence cutoff for the error ellipses
#' @param show.numbers logical flag (\code{TRUE} to show grain numbers)
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
    fit <- yorkfit(x$X,x$sX,x$Y,x$sY,x$rXY)
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
#' @param exterr propagate external sources of uncertainty (J, decay constant)?
#' @return
#' if \code{plot=FALSE}, returns a list with the following items:
#' \describe{
#' \item{a}{the intercept of the straight line fit and its standard error} 
#' \item{b}{the slope of the fit and its standard error}
#' \item{y0}{this either equals \code{a} or, if \code{x} has class
#' \code{ArAr}, the atmospheric \eqn{^{40}}Ar/\eqn{^{36}}Ar ratio and
#' its standard error}
#' \item{age}{the \eqn{^{40}}Ar/\eqn{^{39}}Ar, Re-Os, Rb-Sr or Sm-Nd
#' age and its standard error}
#' }
#' @examples
#' data(examples)
#' isochron(examples$ArAr)
#' @rdname isochron
#' @export
isochron.ArAr <- function(x,xlim=NA,ylim=NA,alpha=0.05,sigdig=2,
                          show.numbers=FALSE,ellipse.col=rgb(0,1,0,0.5),
                          inverse=TRUE,line.col='red',lwd=2,plot=TRUE,
                          exterr=TRUE,...){
    d <- data2york(x,get.selection(x,inverse))
    if (inverse){
        fit <- yorkfit(d$Y,d$sY,d$X,d$sX,-d$rXY) # X and Y reversed!
        x0 <- 1/fit$a[1]
        sx0 <- fit$a[2]/fit$a[1]^2
        y0 <- -fit$b[1]/fit$a[1]
        sy0 <- y0*sqrt((fit$a[2]/fit$a[1])^2 + (fit$b[2]/fit$b[1])^2)
        tt <- get.ArAr.age(x0,sx0,x$J[1],x$J[2],exterr=exterr)
        x.lab <- expression(paste(""^"39","Ar/"^"40","Ar"))
        y.lab <- expression(paste(""^"36","Ar/"^"40","Ar"))
    } else {
        fit <- yorkfit(d$X,d$sX,d$Y,d$sY,d$rXY)
        y0 <- fit$a[1]
        sy0 <- fit$a[2]
        tt <- get.ArAr.age(fit$b[1],fit$b[2],x$J[1],x$J[2],exterr=exterr)
        x.lab <- expression(paste(""^"39","Ar/"^"36","Ar"))
        y.lab <- expression(paste(""^"40","Ar/"^"36","Ar"))
    }
    out <- fit
    class(out) <- "isochron"
    out$y0 <- c(y0,sy0)
    out$age <- tt
    if (plot) {
        isochron.default(d,xlim=xlim,ylim=ylim,alpha=alpha,
                         show.numbers=show.numbers,
                         ellipse.col=ellipse.col,a=fit$a[1],
                         b=fit$b[1], line.col=line.col,lwd=lwd,
                         title=FALSE)
        tt <- roundit(out$age[1],out$age[2])
        title(isochron.title(out,sigdig=sigdig),xlab=x.lab,ylab=y.lab)
    } else {
        return(out)
    }
}
#' @rdname isochron
#' @export
isochron.RbSr <- function(x,xlim=NA,ylim=NA,alpha=0.05,sigdig=2,
                          show.numbers=FALSE,ellipse.col=rgb(0,1,0,0.5),
                          line.col='red',lwd=2,plot=TRUE,exterr=TRUE,...){
    isochron.PD(x,'Rb87',xlim=xlim, ylim=ylim,alpha=alpha,
                sigdig=sigdig, show.numbers=show.numbers,
                ellipse.col=ellipse.col,line.col=line.col, lwd=lwd,
                plot=plot,exterr=exterr,...)
}
#' @rdname isochron
#' @export
isochron.ReOs <- function(x,xlim=NA,ylim=NA,alpha=0.05,sigdig=2,
                          show.numbers=FALSE,ellipse.col=rgb(0,1,0,0.5),
                          line.col='red',lwd=2,plot=TRUE,exterr=TRUE,...){
    isochron.PD(x,'Re187',xlim=xlim, ylim=ylim,alpha=alpha,
                sigdig=sigdig, show.numbers=show.numbers,
                ellipse.col=ellipse.col,line.col=line.col, lwd=lwd,
                plot=plot,exterr=exterr,...)
}
#' @rdname isochron
#' @export
isochron.SmNd <- function(x,xlim=NA,ylim=NA,alpha=0.05,sigdig=2,
                          show.numbers=FALSE,ellipse.col=rgb(0,1,0,0.5),
                          line.col='red',lwd=2,plot=TRUE,exterr=TRUE,...){
    isochron.PD(x,'Sm147',xlim=xlim,ylim=ylim,alpha=alpha,
                sigdig=sigdig, show.numbers=show.numbers,
                ellipse.col=ellipse.col,line.col=line.col, lwd=lwd,
                plot=plot,exterr=exterr,...)
}
isochron.PD <- function(x,nuclide,xlim=NA,ylim=NA, alpha=0.05,
                        sigdig=2,show.numbers=FALSE,
                        ellipse.col=rgb(0,1,0,0.5),line.col='red',
                        lwd=2,plot=TRUE,exterr=TRUE,...){
    if (identical(nuclide,'Sm147')){
        x.lab <- expression(paste(""^"147","Sm/"^"144","Nd"))
        y.lab <- expression(paste(""^"143","Nd/"^"144","Nd"))
    } else if (identical(nuclide,'Re187')){
        x.lab <- expression(paste(""^"187","Re/"^"188","Os"))
        y.lab <- expression(paste(""^"187","Os/"^"188","Os"))
    } else if (identical(nuclide,'Rb87')){
        x.lab <- expression(paste(""^"87","Rb/"^"86","Sr"))
        y.lab <- expression(paste(""^"87","Sr/"^"86","Sr"))
    }
    X <- ppm2ratios(x,exterr=exterr,common=FALSE)
    fit <- yorkfit(X)
    out <- fit
    class(out) <- "isochron"
    out$y0 <- c(fit$a[1],fit$a[2])
    out$age <- get.PD.age(fit$b[1],fit$b[2],nuclide,exterr=exterr)
    if (plot){
        isochron.default(X,xlim=xlim,ylim=ylim,alpha=alpha,
                         show.numbers=show.numbers,
                         ellipse.col=ellipse.col,a=fit$a[1],
                         b=fit$b[1], line.col=line.col,lwd=lwd,
                         title=FALSE)
        tt <- roundit(out$age[1],out$age[2])
        title(isochron.title(out,sigdig=sigdig),xlab=x.lab,ylab=y.lab)
    } else {
        return(out)
    }
}

get.limits <- function(X,sX){
    minx <- min(X-3*sX,na.rm=TRUE)
    maxx <- max(X+3*sX,na.rm=TRUE)    
    c(minx,maxx)
}

isochron.title <- function(fit,sigdig=2){
    rounded.age <- roundit(fit$age[1],fit$age[2],sigdig=sigdig)
    rounded.intercept <- roundit(fit$y0[1],fit$y0[2],sigdig=sigdig)
    line1 <- substitute('age ='~a%+-%b~'(1'~sigma~'), intercept ='~c%+-%d~'(1'~sigma~')',
                        list(a=rounded.age$x, b=rounded.age$err,
                             c=rounded.intercept$x, d=rounded.intercept$err))
    line2 <- substitute('MSWD ='~a~', p('~chi^2*')='~b,
                        list(a=signif(fit$mswd,sigdig), b=signif(fit$p.value,sigdig)))
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
                        list(a=signif(fit$mswd,sigdig), b=signif(fit$p.value,sigdig)))
    graphics::mtext(line1,line=1)
    graphics::mtext(line2,line=0)
}
