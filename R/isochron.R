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
#' @param alpha confidence cutoff for the error ellipses and
#'     confidence intervals
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
#' @param ... optional arguments to be passed on to the generic plot
#'     function if \code{model=2}
#' @references Nicolaysen, L.O., 1961. Graphic interpretation of
#'     discordant age measurements on metamorphic rocks. Annals of the
#'     New York Academy of Sciences, 91(1), pp.198-206.
#'
#' Ludwig, K.R. and Titterington, D.M., 1994. Calculation of
#'     \eqn{^{230}}Th/U isochrons, ages, and errors. Geochimica et
#'     Cosmochimica Acta, 58(22), pp.5031-5042.
#' @rdname isochron
#' @export
isochron <- function(x,...){ UseMethod("isochron",x) }
#' @rdname isochron
#' @export
isochron.default <- function(x,xlim=NA,ylim=NA,alpha=0.05,sigdig=2,
                             show.numbers=FALSE,levels=NA,
                             ellipse.col=c("#00FF0080","#FF000080"),
                             line.col='red',lwd=2,title=TRUE,model=1,...){
    X <- x[,1:5]
    colnames(X) <- c('X','sX','Y','sY','rXY')
    fit <- regression(X,model=model)
    out <- regression_init(fit,alpha=alpha)
    scatterplot(X,xlim=xlim,ylim=ylim,alpha=alpha,
                show.ellipses=1*(model==1),show.numbers=show.numbers,
                levels=levels,ellipse.col=ellipse.col,a=fit$a[1],
                b=fit$b[1],line.col=line.col,lwd=lwd)
    if (title)
        graphics::title(isochrontitle(out,sigdig=sigdig),xlab='X',ylab='Y')
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
#' @param exterr propagate external sources of uncertainty
#' (J, decay constant)?
#'
#' @return if \code{x} has class \code{PbPb}, \code{ArAr},
#'     \code{RbSr}, \code{SmNd}, \code{ReOs} or \code{LuHf}, or
#'     \code{UThHe}, returns a list with the following items:
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
#' \item{df}{the degrees of freedom of the linear fit (\eqn{df=n-2})}
#'
#' \item{y0}{a four-element list containing:
#'
#' \code{y}: the atmospheric \eqn{^{40}}Ar/\eqn{^{36}}Ar or initial
#' \eqn{^{207}}Pb/\eqn{^{204}}Pb, \eqn{^{187}}Os/\eqn{^{188}}Os,
#' \eqn{^{87}}Sr/\eqn{^{86}}Sr, \eqn{^{143}}Nd/\eqn{^{144}}Nd or
#' \eqn{^{176}}Hf/\eqn{^{177}}Hf ratio.
#'
#' \code{s[y]}: the propagated uncertainty of \code{y}
#'
#' \code{ci[y]}: the \eqn{100(1-\alpha/2)\%} confidence interval for
#' \code{y} given the appropriate degrees of freedom.
#'
#' \code{disp[y]}: the \eqn{100(1-\alpha/2)\%} confidence interval for
#' \code{y} enhanced by \eqn{\sqrt{mswd}} (only applicable if \code{
#' model=1}).
#' }
#' 
#' \item{age}{a four-element list containing:
#'
#' \code{t}: the \eqn{^{207}}Pb/\eqn{^{206}}Pb,
#' \eqn{^{40}}Ar/\eqn{^{39}}Ar, \eqn{^{187}}Os/\eqn{^{187}}Re,
#' \eqn{^{87}}Sr/\eqn{^{86}}Sr, \eqn{^{143}}Nd/\eqn{^{144}}Nd or
#' \eqn{^{176}}Hf/\eqn{^{177}}Hf age.
#'
#' \code{s[t]}: the propagated uncertainty of \code{t}
#'
#' \code{ci[t]}: the \eqn{100(1-\alpha/2)\%} confidence interval for
#' \code{t} given the appropriate degrees of freedom.
#'
#' \code{disp[t]}: the \eqn{100(1-\alpha/2)\%} confidence interval for
#' \code{t} enhanced by \eqn{\sqrt{mswd}} (only applicable if \code{
#' model=1}).  }
#'
#' \item{tfact}{the \eqn{100(1-\alpha/2)\%} percentile of a
#' t-distribution with \code{df} degrees of freedom.}
#'
#' }
#'
#' additionally, if \code{model=1}:
#'
#' \describe{
#'
#' \item{mswd}{the mean square of the residuals (a.k.a
#'     `reduced Chi-square') statistic (omitted if \code{model=2}).}
#'
#' \item{p.value}{the p-value of a Chi-square test for linearity}
#'
#' }
#'
#' OR, if \code{x} has class \code{ThU}:
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
#' \item{df}{the degrees of freedom for the linear fit, i.e. \eqn{(3n-3)} if
#' \code{x$format=1} or \code{x$format=2}, and \eqn{(2n-2)} if
#' \code{x$format=3} or \code{x$format=4}}
#'
#' \item{a}{if \code{type=1}: the \eqn{^{230}}Th/\eqn{^{232}}Th
#' intercept; if \code{type=2}: the \eqn{^{230}}Th/\eqn{^{238}}U
#' intercept; if \code{type=3}: the \eqn{^{234}}Th/\eqn{^{232}}Th
#' intercept; if \code{type=4}: the \eqn{^{234}}Th/\eqn{^{238}}U
#' intercept and its propagated uncertainty.}
#' 
#' \item{b}{if \code{type=1}: the \eqn{^{230}}Th/\eqn{^{238}}U slope;
#' if \code{type=2}: the \eqn{^{230}}Th/\eqn{^{232}}Th slope; if
#' \code{type=3}: the \eqn{^{234}}U/\eqn{^{238}}U slope; if
#' \code{type=4}: the \eqn{^{234}}U/\eqn{^{232}}Th slope and its
#' propagated uncertainty.}
#'
#' \item{cov.ab}{the covariance between \code{a} and \code{b}.}
#' 
#' \item{mswd}{the mean square of the residuals (a.k.a `reduced
#'     Chi-square') statistic.}
#'
#' \item{p.value}{the p-value of a Chi-square test for linearity.}
#'
#' \item{y0}{a four-element vector containing:
#'
#' \code{y}: the initial \eqn{^{234}}U/\eqn{^{238}}U-ratio
#'
#' \code{s[y]}: the propagated uncertainty of \code{y}
#'
#' \code{ci[y]}: the \eqn{100(1-\alpha/2)\%} confidence interval for
#' \code{y}
#'
#' \code{disp[y]}: the \eqn{100(1-\alpha/2)\%} confidence interval for
#' \code{y} enhanced by \eqn{\sqrt{mswd}}.}
#'
#' \item{age}{a four-element vector containing:
#'
#' \code{t}: the initial \eqn{^{234}}U/\eqn{^{238}}U-ratio
#'
#' \code{s[t]}: the propagated uncertainty of \code{t}
#'
#' \code{ci[t]}: the \eqn{100(1-\alpha/2)\%} confidence interval for
#' \code{t}
#'
#' \code{disp[t]}: the \eqn{100(1-\alpha/2)\%} confidence interval for
#' \code{t} enhanced by \eqn{\sqrt{mswd}}.}
#'
#' \item{d}{a matrix with the following columns: the X-variable for
#' the isochron plot, the analytical uncertainty of X, the Y-variable
#' for the isochron plot, the analytical uncertainty of Y, and the
#' correlation coefficient between X and Y.}
#'
#' \item{xlab}{the x-label of the isochron plot}
#'
#' \item{ylab}{the y-label of the isochron plot}
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
    out <- isochron_init(fit,alpha=alpha)
    a <- fit$a['a']
    sa <- fit$a['s[a]']
    b <- fit$b['b']
    sb <- fit$b['s[b]']
    if (inverse) {
        R09 <- -b/a
        sR09 <- R09*sqrt((sa/a)^2+(sb/b)^2-2*sa*sb*fit$cov.ab)
        out$y0['y'] <- 1/a
        out$y0['s[y]'] <- sa/a^2
        x.lab <- expression(paste(""^"39","Ar/"^"40","Ar"))
        y.lab <- expression(paste(""^"36","Ar/"^"40","Ar"))
    } else {
        R09 <- b
        sR09 <- sb
        out$y0['y'] <- a
        out$y0['s[y]'] <- sa
        x.lab <- expression(paste(""^"39","Ar/"^"36","Ar"))
        y.lab <- expression(paste(""^"40","Ar/"^"36","Ar"))
    }
    out$age[c('t','s[t]')] <-
        get.ArAr.age(R09,sR09,x$J[1],x$J[2],exterr=exterr)
    out$y0['ci[y]'] <- out$tfact*out$y0['s[y]']
    out$age['ci[t]'] <- out$tfact*out$age['s[t]']
    if (model==1){
        out$y0['disp[y]'] <- sqrt(out$mswd)*out$y0['ci[y]']
        out$age['disp[t]'] <-
            out$tfact*get.ArAr.age(R09,sqrt(out$mswd)*sR09,
                                   x$J[1],x$J[2],exterr=exterr)[2]
    }
    show.ellipses <- (model!=2)
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
    out <- isochron_init(fit,alpha=alpha)
    if (inverse){
        R76 <- fit$a
        out$y0[c('y','s[y]')] <- fit$b
        x.lab <- expression(paste(""^"204","Pb/"^"206","Pb"))
        y.lab <- expression(paste(""^"207","Pb/"^"206","Pb"))
    } else {
        R76 <- fit$b
        out$y0[c('y','s[y]')] <- fit$a
        x.lab <- expression(paste(""^"206","Pb/"^"204","Pb"))
        y.lab <- expression(paste(""^"207","Pb/"^"204","Pb"))
    }
    out$age[c('t','s[t]')] <-
        get.Pb207Pb206.age(R76[1],R76[2],exterr=exterr)
    out$y0['ci[y]'] <- out$tfact*out$y0['s[y]']
    out$age['ci[t]'] <- out$tfact*out$age['s[t]']
    if (model==1){
        out$y0['disp[y]'] <- sqrt(out$mswd)*out$y0['ci[y]']
        out$age['disp[t]'] <-
            out$tfact*get.Pb207Pb206.age(R76[1],sqrt(out$mswd)*R76[2],
                                         exterr=exterr)[2]
    }
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
        out <- isochron_ThU_3D(x,type=type,model=model,
                               exterr=exterr,alpha=alpha)
        intercept.type <- 'Th-U-3D'
    } else if (x$format %in% c(3,4)){
        out <- isochron_ThU_2D(x,type=type,model=model,
                               exterr=exterr,alpha=alpha)
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
    out <- isochron_init(fit,alpha=alpha)
    out$y0[c('y','s[y]')] <- fit$a
    out$age[c('t','s[t]')] <- fit$b
    out$y0['ci[y]'] <- out$tfact*out$y0['s[y]']
    out$age['ci[t]'] <- out$tfact*out$age['s[t]']
    if (model==1){
        out$y0['disp[y]'] <- out$tfact*sqrt(out$mswd)*out$y0['s[y]']
        out$age['disp[t]'] <- out$tfact*sqrt(out$mswd)*out$age['s[t]']
    }
    if (plot) {
        scatterplot(d,xlim=xlim,ylim=ylim,alpha=alpha,
                    show.ellipses=2*(model==1),show.numbers=show.numbers,
                    a=fit$a[1],b=fit$b[1],line.col=line.col,lwd=lwd,...)
        graphics::title(isochrontitle(out,sigdig=sigdig,type='U-Th-He'),
                        xlab="P",ylab="He")
    }
    out
}

isochron_ThU_3D <- function(x,type=2,model=1,
                            exterr=TRUE,alpha=0.05){
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
    out <- isochron_init(fit,alpha=alpha)
    out$a <- c(fit$par[ia],sqrt(fit$cov[ia,ia]))
    out$b <- c(fit$par[ib],sqrt(fit$cov[ib,ib]))
    out$cov.ab <- fit$cov[ia,ib]
    tst <- get.ThU.age(fit$par[i08],sqrt(fit$cov[i08,i08]),
                       fit$par[i48],sqrt(fit$cov[i48,i48]),
                       fit$cov[i48,i08],exterr=exterr)
    out$age['t'] <- tst['t']
    out$y0['y'] <- tst['48_0']
    out$age['s[t]'] <- tst['s[t]']
    out$y0['s[y]'] <- tst['s[48_0]']
    out$age['ci[t]'] <- out$tfact*out$age['s[t]']
    out$y0['ci[y]'] <- out$tfact*out$y0['s[y]']
    if (model==1 && fit$mswd>1){
        tdispt <- get.ThU.age(fit$par[i08],
                              sqrt(out$mswd)*sqrt(fit$cov[i08,i08]),
                              fit$par[i48],
                              sqrt(out$mswd)*sqrt(fit$cov[i48,i48]),
                              out$mswd*fit$cov[i48,i08],
                              exterr=exterr)
        out$age['disp[t]'] <- out$tfact*tdispt['s[t]']
        out$y0['disp[y]'] <- out$tfact*tdispt['s[48_0]']
    }
    out$xlab <- xlab
    out$ylab <- ylab
    out$d <- d[,id]
    out
}
isochron_ThU_2D <- function(x,type=2,model=1,
                            exterr=TRUE,alpha=0.05){
    d <- data2york(x,type=type)
    fit <- regression(d,model=model,type="york")
    if (type==1){
        Th230U238 <- fit$b
        Th230Th232 <- fit$a
        xlab <- expression(paste(""^"238","U/"^"232","Th"))
        ylab <- expression(paste(""^"230","Th/"^"232","Th"))
    } else if (type==2) {
        Th230U238 <- fit$a
        Th230Th232 <- fit$b
        xlab <- expression(paste(""^"232","Th/"^"238","U"))
        ylab <- expression(paste(""^"230","Th/"^"238","U"))
    }
    out <- isochron_init(fit,alpha=alpha)
    out$age[c('t','s[t]')] <-
        get.ThU.age(Th230U238[1],Th230U238[2],
                    exterr=exterr)[c('t','s[t]')]
    out$y0[c('y','s[y]')] <-
        get.Th230Th232_0x(out$age['t'],Th230Th232[1],Th230Th232[2])
    out$age['ci[t]'] <-
        out$tfact*out$age['s[t]']
    out$y0['ci[y]'] <-
        out$tfact*out$y0['s[y]']
    if (model==1 && out$mswd>1){
        out$age['disp[t]'] <-
            out$tfact*get.ThU.age(Th230U238[1],
                                  sqrt(out$mswd)*Th230U238[2],
                                  exterr=exterr)['s[t]']
        out$y0['disp[y]'] <-
            out$tfact*get.Th230Th232_0x(out$age['t'],Th230Th232[1],
                                        sqrt(out$mswd)*Th230Th232[2])[2]
    }
    out$xlab <- xlab
    out$ylab <- ylab
    out$d <- d
    out
}

isochron_PD <- function(x,nuclide,xlim=NA,ylim=NA,alpha=0.05,
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
    out <- isochron_init(fit,alpha=alpha)
    out$y0[c('y','s[y]')] <- fit$a
    out$age[c('t','s[t]')] <- get.PD.age(fit$b['b'],
                   fit$b['s[b]'],nuclide,exterr=exterr)
    out$y0['ci[y]'] <- out$tfact*out$y0['s[y]']
    out$age['ci[t]'] <- out$tfact*out$age['s[t]']
    if (model==1){
        out$y0['disp[y]'] <- out$tfact*sqrt(out$mswd)*out$y0['s[y]']
        out$age['disp[t]'] <- out$tfact*get.PD.age(fit$b['b'],
                              sqrt(out$mswd)*fit$b['s[b]'],
                              nuclide,exterr=exterr)[2]
    }
    if (plot){
        scatterplot(d,xlim=xlim,ylim=ylim,alpha=alpha,
                    show.ellipses=1*(model==1),
                    show.numbers=show.numbers,levels=levels,
                    ellipse.col=ellipse.col, a=fit$a[1],b=fit$b[1],
                    line.col=line.col,lwd=lwd,...)
        graphics::title(isochrontitle(out,sigdig=sigdig,type='PD'),
                        xlab=x.lab,ylab=y.lab)
    }
    out
}

isochron_init <- function(fit,alpha=0.05){
    out <- fit
    out$age <- rep(NA,4)
    out$y0 <- rep(NA,4)
    out$tfact <- stats::qt(1-alpha/2,out$df)
    names(out$age) <- c('t','s[t]','ci[t]','disp[t]')
    names(out$y0) <- c('y','s[y]','ci[y]','disp[y]')
    class(out) <- "isochron"
    out
}
regression_init <- function(fit,alpha=0.05){
    out <- fit
    out$a <- rep(NA,4)
    out$b <- rep(NA,4)
    out$tfact <- stats::qt(1-alpha/2,out$df)
    names(out$a) <- c('a','s[a]','ci[a]','disp[a]')
    names(out$b) <- c('b','s[b]','ci[b]','disp[b]')
    out$a[c('a','s[a]')] <- fit$a[c('a','s[a]')]
    out$b[c('b','s[b]')] <- fit$b[c('b','s[b]')]
    out$a['ci[a]'] <- out$tfact*fit$a['s[a]']
    out$a['disp[a]'] <- out$tfact*sqrt(fit$mswd)*fit$a['s[a]']
    out$b['ci[b]'] <- out$tfact*fit$b['s[b]']
    out$b['disp[b]'] <- out$tfact*sqrt(fit$mswd)*fit$b['s[b]']
    class(out) <- "isochron"
    out
}

get.limits <- function(X,sX){
    minx <- min(X-3*sX,na.rm=TRUE)
    maxx <- max(X+3*sX,na.rm=TRUE)    
    c(minx,maxx)
}

isochrontitle <- function(fit,sigdig=2,type=NA){
    if (fit$model!=2 && fit$mswd>1) args <- quote(a%+-%b~'|'~c~'|'~d)
    else args <- quote(a%+-%b~'|'~c)
    if (is.na(type)){
        intercept <- roundit(fit$a[1],fit$a[2:4],sigdig=sigdig)
        slope <- roundit(fit$b[1],fit$b[2:4],sigdig=sigdig)
        expr1 <- quote('slope =')
        expr2 <- quote('intercept =')
        list1 <- list(a=slope[1],
                      b=slope[2],
                      c=slope[3])
        list2 <- list(a=intercept[1],
                      b=intercept[2],
                      c=intercept[3])
        if (fit$mswd>1){
            list1$d <- slope[4]
            list2$d <- intercept[4]
        }
    } else {
        rounded.age <- roundit(fit$age[1],fit$age[2:4],sigdig=sigdig)
        rounded.intercept <- roundit(fit$y0[1],fit$y0[2:4],sigdig=sigdig)
        expr1 <- quote('age =')
        list1 <- list(a=rounded.age[1],
                      b=rounded.age[2],
                      c=rounded.age[3])
        list2 <- list(a=rounded.intercept[1],
                      b=rounded.intercept[2],
                      c=rounded.intercept[3])
        if (fit$model!=2 && fit$mswd>1){
            list1$d <- rounded.age[4]
            list2$d <- rounded.intercept[4]
        }
        if (identical(type,'Ar-Ar'))
            expr2 <- quote('('^40*'Ar/'^36*'Ar)'[o]~'=')
        else if (identical(type,'Pb-Pb'))
            expr2 <- quote('('^207*'Pb/'^204*'Pb)'[o]~'=')
        else if (identical(type,'Th-U-3D'))
            expr2 <- quote('('^234*'U/'^238*'U)'[o]~'=')
        else if (identical(type,'Th-U-2D'))
            expr2 <- quote('('^230*'Th/'^232*'Th)'[o]^x*~'=')
        else
            expr2 <- quote('intercept =')
    }
    call1 <- substitute(e~a,list(e=expr1,a=args))
    call2 <- substitute(e~a,list(e=expr2,a=args))
    line1 <- do.call(substitute,list(eval(call1),list1))
    line2 <- do.call(substitute,list(eval(call2),list2))
    if (fit$model==1){
        line3 <- substitute('MSWD ='~a~', p('~chi^2*')='~b,
                            list(a=signif(fit$mswd,sigdig),
                                 b=signif(fit$p.value,sigdig)))
        graphics::mtext(line1,line=2)
        graphics::mtext(line2,line=1)
        graphics::mtext(line3,line=0)
    } else {
        graphics::mtext(line1,line=1)
        graphics::mtext(line2,line=0)
    }
}
