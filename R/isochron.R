#' Calculate and plot isochrons
#'
#' Plots cogenetic Ar-Ar, Pb-Pb, Rb-Sr, Sm-Nd, Re-Os, Lu-Hf, U-Th-He or Th-U
#' data as X-Y scatterplots, fits an isochron curve through them using
#' the \code{york} function, and computes the corresponding isochron
#' age, including decay constant uncertainties.
#'
#' @param x EITHER a matrix with the following five columns:
#'
#' \describe{
#'
#' \item{X}{the x-variable}
#'
#' \item{sX}{the standard error of \code{X}}
#'
#' \item{Y}{the y-variable}
#'
#' \item{sY}{the standard error of \code{Y}}
#'
#' \item{rXY}{the correlation coefficient of \code{X} and \code{Y}}
#'
#' }
#'
#' OR
#'
#' an object of class \code{ArAr}, \code{PbPb}, \code{ReOs},
#' \code{RbSr}, \code{SmNd}, \code{LuHf}, \code{UThHe} or \code{ThU}.
#'
#' @param xlim 2-element vector with the plot limits of the x-axis
#'
#' @param ylim 2-element vector with the plot limits of the y-axis
#'
#' @param alpha confidence cutoff for the error ellipses
#'
#' @param show.numbers logical flag (\code{TRUE} to show grain numbers)
#'
#' @param sigdig the number of significant digits of the numerical
#'     values reported in the title of the graphical output
#'
#' @param levels a vector with additional values to be displayed as
#'     different background colours within the error ellipses.
#'
#' @param ellipse.col a vector of two background colours for the error
#'     ellipses. If \code{levels=NA}, then only the first colour will
#'     be used. If \code{levels} is a vector of numbers, then
#'     \code{ellipse.col} is used to construct a colour ramp.
#'
#' @param line.col colour of the isochron line
#'
#' @param lwd line width
#'
#' @param title add a title to the plot?
#'
#' @param model construct the isochron using either:
#'
#' \enumerate{
#'
#' \item{Error weighted least squares regression}
#'
#' \item{Ordinary least squares regression}
#'
#' }
#'
#' @param ... optional arguments to be passed on to the
#' generic plot function if \code{model=2}
#' @rdname isochron
#' @export
isochron <- function(x,...){ UseMethod("isochron",x) }
#' @rdname isochron
#' @export
isochron.default <- function(x,xlim=NA,ylim=NA,alpha=0.05, sigdig=2,
                             show.numbers=FALSE,levels=NA,
                             ellipse.col=c("#00FF0080","#FF000080"),
                             line.col='red',lwd=2,title=TRUE,model=1,...){
    X <- x[,1:5]
    colnames(X) <- c('X','sX','Y','sY','rXY')
    fit <- regression(X,model=model)
    scatterplot(X,xlim=xlim,ylim=ylim,alpha=alpha,
                show.ellipses=1*(model==1),show.numbers=show.numbers,
                levels=levels,ellipse.col=ellipse.col,
                a=fit$a[1],b=fit$b[1],line.col=line.col,lwd=lwd)
    if (title)
        graphics::title(isochrontitle(fit,sigdig=sigdig),xlab='X',ylab='Y')
}
#' @param plot if \code{FALSE}, suppresses the graphical output
#'
#' @param inverse if \code{TRUE} and \code{x} has class \code{ArAr},
#'     plots \eqn{^{36}}Ar/\eqn{^{40}}Ar
#'     vs. \eqn{^{39}}Ar/\eqn{^{40}}Ar.
#'
#' if \code{TRUE} and \code{x} has class \code{PbPb}, plots
#' \eqn{^{207}}Pb/\eqn{^{206}}Pb vs. \eqn{^{204}}Pb/\eqn{^{206}}Pb.
#'
#' @param exterr propagate external sources of uncertainty (J, decay constant)?
#'
#' @return if \code{x} has class \code{PbPb}, \code{ArAr},
#'     \code{RbSr}, \code{SmNd}, \code{ReOs} or \code{LuHf},
#'     \code{ThU}, or \code{UThHe}, returns a list with the following
#'     items:
#'
#' \describe{
#'
#' \item{a}{the intercept of the straight line fit and its standard
#' error.}
#' 
#' \item{b}{the slope of the fit and its standard error.}
#'
#' \item{cov.ab}{the covariance of the slope and intercept}
#' 
#' \item{mswd}{the mean square of the residuals (a.k.a
#'     `reduced Chi-square') statistic}
#'
#' \item{p.value}{the p-value of a Chi-square test for linearity}
#'
#' \item{y0}{the atmospheric \eqn{^{40}}Ar/\eqn{^{36}}Ar or initial
#' \eqn{^{207}}Pb/\eqn{^{204}}Pb, \eqn{^{187}}Os/\eqn{^{188}}Os,
#' \eqn{^{87}}Sr/\eqn{^{86}}Sr, \eqn{^{143}}Nd/\eqn{^{144}}Nd or
#' \eqn{^{176}}Hf/\eqn{^{177}}Hf ratio and its standard error.}
#' 
#' \item{age}{the \eqn{^{207}}Pb/\eqn{^{206}}Pb,
#' \eqn{^{40}}Ar/\eqn{^{39}}Ar, \eqn{^{187}}Os/\eqn{^{187}}Re,
#' \eqn{^{87}}Sr/\eqn{^{86}}Sr, \eqn{^{143}}Nd/\eqn{^{144}}Nd or
#' \eqn{^{176}}Hf/\eqn{^{177}}Hf age and its standard error.}
#'
#' }
#'
#' if \code{plot=FALSE}, and \code{x} has class \code{ThU}:
#'
#' \describe{
#'
#' \item{par}{if \code{type=1} or \code{type=3}: the best fitting
#' \eqn{^{230}}Th/\eqn{^{232}}Th intercept,
#' \eqn{^{230}}Th/\eqn{^{238}}U slope, \eqn{^{234}}U/\eqn{^{232}}Th
#' intercept and \eqn{^{234}}U/\eqn{^{238}}U slope, OR, if
#' \code{type=2} or \code{type=4}: the best fitting
#' \eqn{^{234}}U/\eqn{^{238}}U intercept,
#' \eqn{^{230}}Th/\eqn{^{232}}Th slope, \eqn{^{234}}U/\eqn{^{238}}U
#' intercept and \eqn{^{234}}U/\eqn{^{232}}Th slope.  }
#'
#' \item{cov}{the covariance matrix of \code{par}.}
#'
#' \item{a}{if \code{type=1}: the \eqn{^{230}}Th/\eqn{^{232}}Th
#' intercept; if \code{type=2}: the \eqn{^{230}}Th/\eqn{^{238}}U
#' intercept; if \code{type=3}: the \eqn{^{234}}Th/\eqn{^{232}}Th
#' intercept; if \code{type=4}: the \eqn{^{234}}Th/\eqn{^{238}}U
#' intercept.}
#' 
#' \item{b}{if \code{type=1}: the \eqn{^{230}}Th/\eqn{^{238}}U slope;
#' if \code{type=2}: the \eqn{^{230}}Th/\eqn{^{232}}Th slope; if
#' \code{type=3}: the \eqn{^{234}}U/\eqn{^{238}}U slope; if
#' \code{type=4}: the \eqn{^{234}}U/\eqn{^{232}}Th slope.}
#'
#' \item{cov.ab}{the covariance between \code{a} and \code{b}.}
#' 
#' \item{mswd}{the mean square of the residuals (a.k.a `reduced
#'     Chi-square') statistic.}
#'
#' \item{p.value}{the p-value of a Chi-square test for linearity.}
#'
#' \item{y0}{the initial \eqn{^{234}}U/\eqn{^{238}}U-ratio and its
#' standard error.}
#'
#' \item{age}{the Th-U isochron age and its standard error.}
#'
#' }
#'
#' @examples
#' data(examples)
#' isochron(examples$ArAr)
#'
#' @rdname isochron
#' @export
isochron.ArAr <- function(x,xlim=NA,ylim=NA,alpha=0.05,sigdig=2,
                          show.numbers=FALSE,levels=NA,
                          ellipse.col=c("#00FF0080","#FF000080"),
                          inverse=TRUE,line.col='red',lwd=2,plot=TRUE,
                          exterr=TRUE,model=1,...){
    d <- data2york(x,inverse=inverse)
    fit <- regression(d,model=model)
    if (inverse){
        x0 <- -fit$b[1]/fit$a[1]
        sx0 <- x0*sqrt((fit$a[2]/fit$a[1])^2 + (fit$b[2]/fit$b[1])^2 -
                       2*fit$a[2]*fit$b[2]*fit$cov.ab)
        y0 <- 1/fit$a[1]
        sy0 <- fit$a[2]/fit$a[1]^2
        tt <- get.ArAr.age(x0,sx0,x$J[1],x$J[2],exterr=exterr)
        x.lab <- expression(paste(""^"39","Ar/"^"40","Ar"))
        y.lab <- expression(paste(""^"36","Ar/"^"40","Ar"))
    } else {
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
    show.ellipses <- (model != 2)
    if (plot) {
        scatterplot(d,xlim=xlim,ylim=ylim,alpha=alpha,
                    show.ellipses=show.ellipses,
                    show.numbers=show.numbers,levels=levels,
                    ellipse.col=ellipse.col,a=fit$a[1],b=fit$b[1],
                    line.col=line.col,lwd=lwd,...)
        graphics::title(isochrontitle(out,sigdig=sigdig,type='Ar-Ar'),
                        xlab=x.lab,ylab=y.lab)
    }
    out
}
#' @rdname isochron
#' @export
isochron.PbPb <- function(x,xlim=NA,ylim=NA,alpha=0.05,sigdig=2,
                          show.numbers=FALSE,levels=NA,
                          ellipse.col=c("#00FF0080","#FF000080"),
                          inverse=TRUE,line.col='red',lwd=2,plot=TRUE,
                          exterr=TRUE,model=1,...){
    d <- data2york(x,inverse=inverse)
    fit <- regression(d,model=model)
    if (inverse){
        y0 <- fit$b[1]
        sy0 <- fit$b[2]
        tt <- get.Pb207Pb206.age(fit$a[1],fit$a[2],exterr=exterr)
        x.lab <- expression(paste(""^"204","Pb/"^"206","Pb"))
        y.lab <- expression(paste(""^"207","Pb/"^"206","Pb"))
    } else {
        y0 <- fit$a[1]
        sy0 <- fit$a[2]
        tt <- get.Pb207Pb206.age(fit$b[1],fit$b[2],exterr=exterr)
        x.lab <- expression(paste(""^"206","Pb/"^"204","Pb"))
        y.lab <- expression(paste(""^"207","Pb/"^"204","Pb"))
    }
    out <- fit
    class(out) <- "isochron"
    out$y0 <- c(y0,sy0)
    out$age <- tt
    if (plot) {
        scatterplot(d,xlim=xlim,ylim=ylim,alpha=alpha,
                    show.ellipses=1*(model==1),
                    show.numbers=show.numbers,levels=levels,
                    ellipse.col=ellipse.col,a=fit$a[1],b=fit$b[1],
                    line.col=line.col,lwd=lwd,...)
        graphics::title(isochrontitle(out,sigdig=sigdig,type='Pb-Pb'),
                        xlab=x.lab,ylab=y.lab)
    }
    out
}
#' @rdname isochron
#' @export
isochron.RbSr <- function(x,xlim=NA,ylim=NA,alpha=0.05,sigdig=2,
                          show.numbers=FALSE,levels=NA,
                          ellipse.col=c("#00FF0080","#FF000080"),
                          line.col='red',lwd=2,plot=TRUE,exterr=TRUE,
                          model=1,...){
    isochron_PD(x,'Rb87',xlim=xlim,ylim=ylim,alpha=alpha,
                sigdig=sigdig,show.numbers=show.numbers,
                levels=levels,ellipse.col=ellipse.col,
                line.col=line.col,lwd=lwd,plot=plot,exterr=exterr,
                model=model,...)
}
#' @rdname isochron
#' @export
isochron.ReOs <- function(x,xlim=NA,ylim=NA,alpha=0.05,sigdig=2,
                          show.numbers=FALSE,levels=NA,
                          ellipse.col=c("#00FF0080","#FF000080"),
                          line.col='red',lwd=2,plot=TRUE,exterr=TRUE,
                          model=1,...){
    isochron_PD(x,'Re187',xlim=xlim,ylim=ylim,alpha=alpha,
                sigdig=sigdig,show.numbers=show.numbers,
                levels=levels,ellipse.col=ellipse.col,
                line.col=line.col,lwd=lwd,plot=plot,exterr=exterr,
                model=model,...)
}
#' @rdname isochron
#' @export
isochron.SmNd <- function(x,xlim=NA,ylim=NA,alpha=0.05,sigdig=2,
                          show.numbers=FALSE,levels=NA,
                          ellipse.col=c("#00FF0080","#FF000080"),
                          line.col='red',lwd=2,plot=TRUE,exterr=TRUE,
                          model=1,...){
    isochron_PD(x,'Sm147',xlim=xlim,ylim=ylim,alpha=alpha,
                sigdig=sigdig,show.numbers=show.numbers,
                levels=levels,ellipse.col=ellipse.col,
                line.col=line.col,lwd=lwd,plot=plot,exterr=exterr,
                model=model,...)
}
#' @rdname isochron
#' @export
isochron.LuHf <- function(x,xlim=NA,ylim=NA,alpha=0.05,sigdig=2,
                          show.numbers=FALSE,levels=NA,
                          ellipse.col=c("#00FF0080","#FF000080"),
                          line.col='red',lwd=2,plot=TRUE,exterr=TRUE,
                          model=1,...){
    isochron_PD(x,'Lu176',xlim=xlim,ylim=ylim,alpha=alpha,
                sigdig=sigdig,show.numbers=show.numbers,
                levels=levels,ellipse.col=ellipse.col,
                line.col=line.col,lwd=lwd,plot=plot,exterr=exterr,
                model=model,...)
}
#' @param type following the classification of
#' Ludwig and Titterington (1994), one of either:
#' 
#' \enumerate{
#' 
#' \item `Rosholt type-II' isochron, setting out
#' \eqn{^{230}}Th/\eqn{^{232}}Th vs. \eqn{^{238}}U/\eqn{^{232}}Th
#' 
#' \item `Osmond type-II' isochron, setting out \eqn{^{230}}Th/\eqn{^{238}}U
#' vs. \eqn{^{232}}Th/\eqn{^{238}}U
#'
#' \item `Rosholt type-II' isochron, setting out \eqn{^{234}}U/\eqn{^{232}}Th
#' vs. \eqn{^{238}}U/\eqn{^{232}}Th
#' 
#' \item `Osmond type-II' isochron, setting out \eqn{^{234}}U/\eqn{^{238}}U
#' vs. \eqn{^{232}}Th/\eqn{^{238}}U
#'
#' }
#' @rdname isochron
#' @export
isochron.ThU <- function (x,type=2,xlim=NA,ylim=NA,alpha=0.05,
                          sigdig=2,show.numbers=FALSE,levels=NA,
                          ellipse.col=c("#00FF0080","#FF000080"),
                          line.col='red',lwd=2,plot=TRUE,exterr=TRUE,
                          model=1,...){
    if (x$format %in% c(1,2)){
        out <- isochron_ThU_3D(x,type=type,model=model,exterr=exterr)
        intercept.type <- 'Th-U-3D'
    } else if (x$format %in% c(3,4)){
        out <- isochron_ThU_2D(x,type=type,model=model,exterr=exterr)
        intercept.type <- 'Th-U-2D'
    }
    if (plot){
        scatterplot(out$d,xlim=xlim,ylim=ylim,alpha=alpha,
                    show.ellipses=1*(model==1),
                    show.numbers=show.numbers,levels=levels,
                    ellipse.col=ellipse.col,a=out$a[1],b=out$b[1],
                    line.col=line.col,lwd=lwd,...)
        graphics::title(isochrontitle(out,sigdig=sigdig,type=intercept.type),
                        xlab=out$xlab,ylab=out$ylab)
    }
    out
}
#' @rdname isochron
#' @export
isochron.UThHe <- function(x,xlim=NA,ylim=NA,alpha=0.05,sigdig=2,
                           show.numbers=FALSE,line.col='red',lwd=2,
                           plot=TRUE,model=1,...){
    d <- data2york(x)
    fit <- regression(d,model=model)
    out <- fit
    class(out) <- "isochron"
    out$y0 <- fit$a
    out$age <- fit$b
    if (plot) {
        scatterplot(d,xlim=xlim,ylim=ylim,alpha=alpha,
                    show.ellipses=2*(model==1),show.numbers=show.numbers,
                    a=fit$a[1],b=fit$b[1],line.col=line.col,lwd=lwd,...)
        graphics::title(isochrontitle(out,sigdig=sigdig,type='U-Th-He'),
                        xlab="P",ylab="He")
    }
    out
}

isochron_ThU_3D <- function(x,type=2,model=1,exterr=TRUE){
    if (type == 1){
        osmond <- FALSE
        ia <- 'a'
        ib <- 'b'
        i48 <- 'b'
        i08 <- 'B'
        id <- c('X','sX','Y','sY','rXY')
        xlab <- expression(paste(""^"238","U/"^"232","Th"))
        ylab <- expression(paste(""^"230","Th/"^"232","Th"))
    } else if (type == 2){
        osmond <- TRUE
        ia <- 'A'
        ib <- 'B'
        i48 <- 'a'
        i08 <- 'A'
        id <- c('X','sX','Z','sZ','rXZ')
        xlab <- expression(paste(""^"232","Th/"^"238","U"))
        ylab <- expression(paste(""^"230","Th/"^"238","U"))
    } else if (type == 3){
        osmond <- FALSE
        ia <- 'A'
        ib <- 'B'
        i48 <- 'b'
        i08 <- 'B'
        id <- c('X','sX','Z','sZ','rXZ')
        xlab <- expression(paste(""^"238","U/"^"232","Th"))
        ylab <- expression(paste(""^"234","U/"^"232","Th"))
    } else {
        osmond <- TRUE
        ia <- 'a'
        ib <- 'b'
        i48 <- 'a'
        i08 <- 'A'
        id <- c('X','sX','Y','sY','rXY')
        xlab <- expression(paste(""^"232","Th/"^"238","U"))
        ylab <- expression(paste(""^"234","U/"^"238","U"))
    }
    d <- data2tit(x,osmond=osmond)
    fit <- regression(d,model=model,type="titterington")
    out <- fit
    class(out) <- "isochron"
    out$a <- c(fit$par[ia],sqrt(fit$cov[ia,ia]))
    out$b <- c(fit$par[ib],sqrt(fit$cov[ib,ib]))
    out$cov.ab <- fit$cov[ia,ib]
    tt <- get.ThU.age(fit$par[i08],sqrt(fit$cov[i08,i08]),
                      fit$par[i48],sqrt(fit$cov[i48,i48]),
                      fit$cov[i48,i08],exterr=exterr)
    out$y0 <- tt[c('48_0','s[48_0]')]
    out$age <- tt[c('t','s[t]')]
    out$xlab <- xlab
    out$ylab <- ylab
    out$d <- d[,id]
    out
}
isochron_ThU_2D <- function(x,type=2,model=1,exterr=TRUE){
    d <- data2york(x,type=type)
    fit <- regression(d,model=model,type="york")
    out <- fit
    class(out) <- "isochron"
    if (type==1){
        Th230U238 <- fit$b
        xlab <- expression(paste(""^"238","U/"^"232","Th"))
        ylab <- expression(paste(""^"230","Th/"^"232","Th"))
    } else if (type==2) {
        Th230U238 <- fit$a
        xlab <- expression(paste(""^"232","Th/"^"238","U"))
        ylab <- expression(paste(""^"230","Th/"^"238","U"))
    }
    tt <- get.ThU.age(Th230U238[1],Th230U238[2],exterr=exterr)
    out$age <- tt[c('t','s[t]')]
    if (type==1)
        out$y0 <- get.Th230Th232_0x(tt,fit$a[1],fit$a[2])
    else
        out$y0 <- get.Th230Th232_0x(tt,fit$b[1],fit$b[2])
    out$xlab <- xlab
    out$ylab <- ylab
    out$d <- d
    out    
}

isochron_PD <- function(x,nuclide,xlim=NA,ylim=NA, alpha=0.05,
                        sigdig=2,show.numbers=FALSE,levels=NA,
                        ellipse.col=c("#00FF0080","#FF000080"),
                        line.col='red',lwd=2,plot=TRUE,exterr=TRUE,
                        model=1,...){
    if (identical(nuclide,'Sm147')){
        x.lab <- expression(paste(""^"147","Sm/"^"144","Nd"))
        y.lab <- expression(paste(""^"143","Nd/"^"144","Nd"))
    } else if (identical(nuclide,'Re187')){
        x.lab <- expression(paste(""^"187","Re/"^"188","Os"))
        y.lab <- expression(paste(""^"187","Os/"^"188","Os"))
    } else if (identical(nuclide,'Rb87')){
        x.lab <- expression(paste(""^"87","Rb/"^"86","Sr"))
        y.lab <- expression(paste(""^"87","Sr/"^"86","Sr"))
    } else if (identical(nuclide,'Lu176')){
        x.lab <- expression(paste(""^"176","Lu/"^"177","Hf"))
        y.lab <- expression(paste(""^"176","Hf/"^"177","Hf"))
    }
    d <- data2york(x,exterr=exterr,common=FALSE)
    fit <- regression(d,model=model)
    out <- fit
    class(out) <- "isochron"
    out$y0 <- c(fit$a[1],fit$a[2])
    out$age <- get.PD.age(fit$b[1],fit$b[2],nuclide,exterr=exterr)
    if (plot){
        scatterplot(d,xlim=xlim,ylim=ylim,alpha=alpha,
                    show.ellipses=1*(model==1),
                    show.numbers=show.numbers,levels=levels,
                    ellipse.col=ellipse.col, a=fit$a[1],b=fit$b[1],
                    line.col=line.col,lwd=lwd,...)
        graphics::title(isochrontitle(out,sigdig=sigdig,type='PD'),
                        xlab=x.lab,ylab=y.lab)
    } else {
        return(out)
    }
}

get.limits <- function(X,sX){
    minx <- min(X-3*sX,na.rm=TRUE)
    maxx <- max(X+3*sX,na.rm=TRUE)    
    c(minx,maxx)
}

isochrontitle <- function(fit,sigdig=2,type=NA){
    if (is.na(type)){
        intercept <- roundit(fit$a[1],fit$a[2],sigdig=sigdig)
        slope <- roundit(fit$b[1],fit$b[2],sigdig=sigdig)
        line1 <- substitute('slope ='~a%+-%b~'(1'~sigma~')',
                            list(a=slope[1], b=slope[2]))
        line2 <- substitute('intercept ='~c%+-%d~'(1'~sigma~')',
                            list(c=intercept[1], d=intercept[2]))
    } else {
        rounded.age <- roundit(fit$age[1],fit$age[2],sigdig=sigdig)
        rounded.intercept <- roundit(fit$y0[1],fit$y0[2],sigdig=sigdig)
        line1 <- substitute('age ='~a%+-%b~'(1'~sigma~')',
                            list(a=rounded.age[1], b=rounded.age[2]))
        if (identical(type,'Ar-Ar')) {
            line2 <- substitute('('^40*'Ar/'^39*'Ar)'[o]~'='~c%+-%d~'(1'~sigma~')',
                                list(c=rounded.intercept[1], d=rounded.intercept[2]))
        } else if (identical(type,'Pb-Pb')) {
            line2 <- substitute('('^207*'Pb/'^204*'Pb)'[o]~'='~c%+-%d~'(1'~sigma~')',
                                list(c=rounded.intercept[1], d=rounded.intercept[2]))
        } else if (identical(type,'Th-U-3D')) {
            line2 <- substitute('('^234*'U/'^238*'U)'[o]~'='~c%+-%d~'(1'~sigma~')',
                                list(c=rounded.intercept[1], d=rounded.intercept[2]))
        } else if (identical(type,'Th-U-2D')) {
            line2 <- substitute('('^230*'Th/'^232*'Th)'[o]^x*~'='~c%+-%d~'(1'~sigma~')',
                                list(c=rounded.intercept[1], d=rounded.intercept[2]))
        } else {
            line2 <- substitute('intercept ='~c%+-%d~'(1'~sigma~')',
                                list(c=rounded.intercept[1], d=rounded.intercept[2]))
        }
    }
    if (fit$model==1){
        line3 <- substitute('MSWD ='~a~', p('~chi^2*')='~b,
                            list(a=signif(fit$mswd,sigdig), b=signif(fit$p.value,sigdig)))
        graphics::mtext(line1,line=2)
        graphics::mtext(line2,line=1)
        graphics::mtext(line3,line=0)
    } else {
        graphics::mtext(line1,line=1)
        graphics::mtext(line2,line=0)
    }
}
