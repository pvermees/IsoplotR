#' Calculate and plot isochrons
#'
#' Plots cogenetic Ar-Ar, K-Ca, Pb-Pb, Rb-Sr, Sm-Nd, Re-Os, Lu-Hf,
#' U-Th-He or Th-U data as X-Y scatterplots, fits an isochron curve
#' through them using the \code{york} function, and computes the
#' corresponding isochron age, including decay constant uncertainties.
#'
#' @details
#' Given several aliquots from a single sample, isochrons allow the
#' non-radiogenic component of the daughter nuclide to be quantified
#' and separated from the radiogenic component. In its simplest form,
#' an isochron is obtained by setting out the amount of radiogenic
#' daughter against the amount of radioactive parent, both normalised
#' to a non-radiogenic isotope of the daughter element, and fitting a
#' straight line through these points by least squares regression
#' (Nicolaysen, 1961). The slope and intercept then yield the
#' radiogenic daughter-parent ratio and the non-radiogenic daughter
#' composition, respectively. There are several ways to fit an
#' isochron.  The easiest of these is ordinary least squares
#' regression, which weighs all data points equally. In the presence
#' of quantifiable analytical uncertainty, it is equally
#' straightforward to use the inverse of the y-errors as weights.  It
#' is significantly more difficult to take into account uncertainties
#' in both the x- and the y-variable (York, 1966). \code{IsoplotR}
#' does so for its U-Th-He isochron calculations. The York (1966)
#' method assumes that the analytical uncertainties of the x- and
#' y-variables are independent from each other. This assumption is
#' rarely met in geochronology.  York (1968) addresses this issue with
#' a bivariate error weighted linear least squares algorithm that
#' accounts for covariant errors in both variables. This algorithm was
#' further improved by York et al. (2004) to ensure consistency with
#' the maximum likelihood approach of Titterington and Halliday
#' (1979).
#'
#' \code{IsoplotR} uses the York et al. (2004) algorithm for its
#' Ar-Ar, K-Ca, Pb-Pb, Rb-Sr, Sm-Nd, Re-Os and Lu-Hf isochrons. The
#' maximum likelihood algorithm of Titterington and Halliday (1979)
#' was generalised from two to three dimensions by Ludwig and
#' Titterington (1994) for U-series disequilibrium dating.  Also this
#' algorithm is implemented in \code{IsoplotR}.  The extent to which
#' the observed scatter in the data can be explained by the analytical
#' uncertainties can be assessed using the Mean Square of the Weighted
#' Deviates (MSWD, McIntyre et al., 1966), which is defined as:
#'
#' \eqn{MSWD = ([X - \hat{X}] \Sigma_{X}^{-1} [X - \hat{X}]^T)/df}
#'
#' where \eqn{X} are the data, \eqn{\hat{X}} are the fitted values,
#' and \eqn{\Sigma_X} is the covariance matrix of \eqn{X}, and \eqn{df
#' = k(n-1)} are the degrees of freedom, where \eqn{k} is the
#' dimensionality of the linear fit. MSWD values that are far smaller
#' or greater than 1 indicate under- or overdispersed measurements,
#' respectively. Underdispersion can be attributed to overestimated
#' analytical uncertainties. \code{IsoplotR} provides three
#' alternative strategies to deal with overdispersed data:
#'
#' \enumerate{
#'
#' \item Attribute the overdispersion to an underestimation of the
#' analytical uncertainties. In this case, the excess scatter can be
#' accounted for by inflating those uncertainties by a \emph{factor}
#' \eqn{\sqrt{MSWD}}.
#'
#' \item Ignore the analytical uncertainties and perform an ordinary
#' least squares regression.
#'
#' \item Attribute the overdispersion to the presence of `geological
#' scatter'.  In this case, the excess scatter can be accounted for by
#' adding an overdispersion \emph{term} that lowers the MSWD to unity.
#'
#' }
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
#' an object of class \code{ArAr}, \code{KCa}, \code{PbPb},
#' \code{ReOs}, \code{RbSr}, \code{SmNd}, \code{LuHf}, \code{UThHe} or
#' \code{ThU}.
#'
#' @param xlim 2-element vector with the x-axis limits
#'
#' @param ylim 2-element vector with the y-axis limits
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
#' @param clabel label for the colour scale
#'
#' @param ellipse.col a vector of two background colours for the error
#'     ellipses. If \code{levels=NA}, then only the first colour will
#'     be used. If \code{levels} is a vector of numbers, then
#'     \code{ellipse.col} is used to construct a colour ramp.
#'
#' @param ci.col the fill colour for the confidence interval of the
#'     intercept and slope.
#'
#' @param line.col colour of the isochron line
#'
#' @param lwd line width
#'
#' @param plot if \code{FALSE}, suppresses the graphical output
#'
#' @param title add a title to the plot?
#'
#' @param model construct the isochron using either:
#'
#' \enumerate{
#'
#' \item{Error-weighted least squares regression}
#'
#' \item{Ordinary least squares regression}
#'
#' \item{Error-weighted least squares with overdispersion term}
#'
#' }
#'
#' @param show.ellipses show the data as:
#'
#' \enumerate{
#'
#' \item points
#'
#' \item error ellipses
#'
#' \item error crosses
#'
#' }
#'
#' @param xlab text label for the horizontal plot axis
#' 
#' @param hide vector with indices of aliquots that should be removed
#'     from the plot.
#' 
#' @param omit vector with indices of aliquots that should be plotted
#'     but omitted from the isochron age calculation.
#' 
#' @param omit.col colour that should be used for the omitted
#'     aliquots.
#' 
#' @param ylab text label for the vertical plot axis
#' 
#' @seealso
#' \code{\link{york}},
#' \code{\link{titterington}},
#' \code{\link{ludwig}}
#' @param ... optional arguments to be passed on to the generic plot
#'     function if \code{model=2}
#' @references
#' Ludwig, K.R. and Titterington, D.M., 1994. Calculation of
#' \eqn{^{230}}Th/U isochrons, ages, and errors. Geochimica et
#' Cosmochimica Acta, 58(22), pp.5031-5042.
#'
#' Nicolaysen, L.O., 1961. Graphic interpretation of discordant age
#' measurements on metamorphic rocks. Annals of the New York Academy
#' of Sciences, 91(1), pp.198-206.
#'
#' Titterington, D.M. and Halliday, A.N., 1979. On the fitting of
#' parallel isochrons and the method of maximum likelihood. Chemical
#' Geology, 26(3), pp.183-195.
#'
#' York, D., 1966. Least-squares fitting of a straight line. Canadian
#' Journal of Physics, 44(5), pp.1079-1086.
#'
#' York, D., 1968. Least squares fitting of a straight line with
#' correlated errors. Earth and Planetary Science Letters, 5,
#' pp.320-324.
#'
#' York, D., Evensen, N.M., Martinez, M.L. and De Basebe Delgado, J., 2004.
#' Unified equations for the slope, intercept, and standard
#' errors of the best straight line. American Journal of Physics,
#' 72(3), pp.367-375.
#'
#' @rdname isochron
#' @export
isochron <- function(x,...){ UseMethod("isochron",x) }
#' @rdname isochron
#' @export
isochron.default <- function(x,xlim=NA,ylim=NA,alpha=0.05,sigdig=2,
                             show.numbers=FALSE,levels=NA,clabel="",
                             ellipse.col=c("#00FF0080","#FF000080"),
                             ci.col='gray80',line.col='black',lwd=1,
                             plot=TRUE,title=TRUE,model=1,
                             show.ellipses=1*(model!=2),xlab='x',
                             ylab='y',hide=NULL,omit=NULL,
                             omit.col=NA,...){
    y <- data2york(x)
    d2calc <- clear(x,hide,omit)
    fit <- regression(d2calc,model=model)
    fit <- regression_init(fit,alpha=alpha)
    fit <- ci_isochron(fit)
    if (plot){
        scatterplot(y,xlim=xlim,ylim=ylim,alpha=alpha,
                    show.ellipses=show.ellipses,
                    show.numbers=show.numbers,levels=levels,
                    clabel=clabel,ellipse.col=ellipse.col,fit=fit,
                    ci.col=ci.col,line.col=line.col,lwd=lwd,
                    hide=hide,omit=omit,omit.col=omit.col,...)
        if (title) graphics::title(isochrontitle(fit,sigdig=sigdig),
                                   xlab=xlab,ylab=ylab)
    }
    invisible(fit)
}
#' @param inverse
#' if \code{FALSE} and \code{x} has class \code{ArAr}, plots
#'     \eqn{^{40}}Ar/\eqn{^{36}}Ar vs. \eqn{^{39}}Ar/\eqn{^{36}}Ar.
#'
#' if \code{FALSE} and \code{x} has class \code{PbPb}, plots
#' \eqn{^{207}}Pb/\eqn{^{204}}Pb vs. \eqn{^{206}}Pb/\eqn{^{204}}Pb.
#'
#' if \code{TRUE} and \code{x} has class \code{ArAr}, plots
#'     \eqn{^{36}}Ar/\eqn{^{40}}Ar vs. \eqn{^{39}}Ar/\eqn{^{40}}Ar.
#'
#' if \code{TRUE} and \code{x} has class \code{PbPb}, plots
#' \eqn{^{207}}Pb/\eqn{^{206}}Pb vs. \eqn{^{204}}Pb/\eqn{^{206}}Pb.
#'
#' @param exterr propagate external sources of uncertainty
#' (J, decay constant)?
#'
#' @return If \code{x} has class \code{PbPb}, \code{ArAr}, \code{KCa},
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
#' \eqn{^{40}}Ca/\eqn{^{44}}Ca, \eqn{^{187}}Os/\eqn{^{188}}Os,
#' \eqn{^{87}}Sr/\eqn{^{87}}Rb, \eqn{^{143}}Nd/\eqn{^{144}}Nd or
#' \eqn{^{176}}Hf/\eqn{^{177}}Hf ratio.
#'
#' \code{s[y]}: the propagated uncertainty of \code{y}
#'
#' \code{ci[y]}: the studentised \eqn{100(1-\alpha)\%} confidence
#' interval for \code{y}.
#'
#' \code{disp[y]}: the studentised \eqn{100(1-\alpha)\%} confidence
#' interval for \code{y} enhanced by \eqn{\sqrt{mswd}} (only
#' applicable if \code{ model=1}).  }
#'
#' \item{age}{a four-element list containing:
#'
#' \code{t}: the \eqn{^{207}}Pb/\eqn{^{206}}Pb,
#' \eqn{^{40}}Ar/\eqn{^{39}}Ar, \eqn{^{40}}K/\eqn{^{40}}Ca,
#' \eqn{^{187}}Os/\eqn{^{187}}Re, \eqn{^{87}}Sr/\eqn{^{87}}Rb,
#' \eqn{^{143}}Nd/\eqn{^{144}}Nd or \eqn{^{176}}Hf/\eqn{^{177}}Hf age.
#'
#' \code{s[t]}: the propagated uncertainty of \code{t}
#'
#' \code{ci[t]}: the studentised \eqn{100(1-\alpha)\%} confidence
#' interval for \code{t}.
#'
#' \code{disp[t]}: the studentised \eqn{100(1-\alpha)\%} confidence
#' interval for \code{t} enhanced by \eqn{\sqrt{mswd}} (only
#' applicable if \code{ model=1}).  }
#'
#' \item{mswd}{the mean square of the residuals (a.k.a `reduced
#'     Chi-square') statistic (omitted if \code{model=2}).}
#'
#' \item{p.value}{the p-value of a Chi-square test for linearity
#' (omitted if \code{model=2})}
#'
#' \item{w}{the overdispersion term, i.e. a three-element vector with
#' the standard deviation of the (assumedly) Normally distributed
#' geological scatter that underlies the measurements, and the lower
#' and upper half-widths of its \eqn{100(1-\alpha)\%} confidence
#' interval (only returned if \code{model=3}).}
#'
#' }
#'
#' OR, if \code{x} has class \code{ThU}:
#'
#' \describe{
#'
#' \item{par}{if \code{x$type=1} or \code{x$type=3}: the best fitting
#' \eqn{^{230}}Th/\eqn{^{232}}Th intercept,
#' \eqn{^{230}}Th/\eqn{^{238}}U slope, \eqn{^{234}}U/\eqn{^{232}}Th
#' intercept and \eqn{^{234}}U/\eqn{^{238}}U slope, OR, if
#' \code{x$type=2} or \code{x$type=4}: the best fitting
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
#' \item{tfact}{the \eqn{100(1-\alpha/2)\%} percentile of a
#' t-distribution with \code{df} degrees of freedom.}
#'
#' \item{y0}{a four-element vector containing:
#'
#' \code{y}: the initial \eqn{^{234}}U/\eqn{^{238}}U-ratio
#'
#' \code{s[y]}: the propagated uncertainty of \code{y}
#'
#' \code{ci[y]}: the studentised \eqn{100(1-\alpha)\%} confidence
#' interval for \code{y}.
#'
#' \code{disp[y]}: the studentised \eqn{100(1-\alpha)\%} confidence
#' interval for \code{y} enhanced by \eqn{\sqrt{mswd}}.}
#'
#' \item{age}{a three (or four) element vector containing:
#'
#' \code{t}: the initial \eqn{^{234}}U/\eqn{^{238}}U-ratio
#'
#' \code{s[t]}: the propagated uncertainty of \code{t}
#'
#' \code{ci[t]}: the studentised \eqn{100(1-\alpha)\%} confidence
#' interval for \code{t}
#'
#' \code{disp[t]}: the studentised \eqn{100(1-\alpha)\%} confidence
#' interval for \code{t} enhanced by \eqn{\sqrt{mswd}} (only reported
#' if \code{model=1}).}
#'
#' \item{w}{the overdispersion term, i.e. a three-element vector with
#' the standard deviation of the (assumedly) Normally distributed
#' geological scatter that underlies the measurements, and the lower
#' and upper half-width of its \eqn{100(1-\alpha)\%} confidence
#' interval (only returned if \code{model=3}).}
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
#' isochron(examples$RbSr)
#'
#' fit <- isochron(examples$ArAr,inverse=FALSE,plot=FALSE)
#'
#' dev.new()
#' isochron(examples$ThU,type=4)
#' @rdname isochron
#' @export
isochron.ArAr <- function(x,xlim=NA,ylim=NA,alpha=0.05,sigdig=2,
                          show.numbers=FALSE,levels=NA,clabel="",
                          ellipse.col=c("#00FF0080","#FF000080"),
                          inverse=TRUE,ci.col='gray80',line.col='black',
                          lwd=1,plot=TRUE,exterr=TRUE,model=1,
                          show.ellipses=1*(model!=2),
                          hide=NULL,omit=NULL,omit.col=NA,...){
    y <- data2york(x,inverse=inverse)
    d2calc <- clear(y,hide,omit)
    fit <- regression(d2calc,model=model)
    out <- isochron_init(fit,alpha=alpha)
    a <- out$a['a']
    sa <- out$a['s[a]']
    b <- out$b['b']
    sb <- out$b['s[b]']
    if (inverse) {
        R09 <- -b/a
        sR09 <- R09*sqrt((sa/a)^2+(sb/b)^2-2*sa*sb*out$cov.ab)
        out$y0['y'] <- 1/a
        out$y0['s[y]'] <- sa/a^2
        x.lab <- quote(''^39*'Ar/'^40*'Ar')
        y.lab <- quote(''^36*'Ar/'^40*'Ar')
    } else {
        R09 <- b
        sR09 <- sb
        out$y0['y'] <- a
        out$y0['s[y]'] <- sa
        x.lab <- quote(''^39*'Ar/'^36*'Ar')
        y.lab <- quote(''^40*'Ar/'^36*'Ar')
    }
    out$displabel <-
        substitute(a*b*c,list(a='(',b=y.lab,c=')-dispersion = '))
    out$y0label <- quote('('^40*'Ar/'^36*'Ar)'[o]*' = ')
    out$age[c('t','s[t]')] <-
        get.ArAr.age(R09,sR09,x$J[1],x$J[2],exterr=exterr)
    out <- ci_isochron(out)
    if (model==1){
        out$age['disp[t]'] <-
            out$fact*get.ArAr.age(R09,sqrt(out$mswd)*sR09,
                                  x$J[1],x$J[2],exterr=exterr)[2]
    }
    if (plot) {
        scatterplot(y,xlim=xlim,ylim=ylim,alpha=alpha,
                    show.ellipses=show.ellipses,
                    show.numbers=show.numbers,levels=levels,
                    clabel=clabel,ellipse.col=ellipse.col,fit=out,
                    ci.col=ci.col,line.col=line.col,lwd=lwd,
                    hide=hide,omit=omit,omit.col=omit.col,...)
        graphics::title(isochrontitle(out,sigdig=sigdig,type='Ar-Ar'),
                        xlab=x.lab,ylab=y.lab)
    }
    invisible(out)
}
#' @rdname isochron
#' @export
isochron.KCa <- function(x,xlim=NA,ylim=NA,alpha=0.05,sigdig=2,
                         show.numbers=FALSE,levels=NA,clabel="",
                         ellipse.col=c("#00FF0080","#FF000080"),
                         ci.col='gray80',line.col='black',
                         lwd=1,plot=TRUE,exterr=TRUE,model=1,
                         show.ellipses=1*(model!=2),
                         hide=NULL,omit=NULL,omit.col=NA,...){
    isochron_PD(x,'K40',xlim=xlim,ylim=ylim,alpha=alpha,
                sigdig=sigdig,show.numbers=show.numbers,
                levels=levels,clabel=clabel,ellipse.col=ellipse.col,
                ci.col=ci.col,line.col=line.col,lwd=lwd,
                plot=plot,exterr=exterr,model=model,
                show.ellipses=show.ellipses,bratio=0.895,
                hide=hide,omit=omit,omit.col=omit.col,...)
}
#' @param growth add Stacey-Kramers Pb-evolution curve to the plot?
#' @rdname isochron
#' @export
isochron.PbPb <- function(x,xlim=NA,ylim=NA,alpha=0.05,sigdig=2,
                          show.numbers=FALSE,levels=NA,clabel="",
                          ellipse.col=c("#00FF0080","#FF000080"),
                          inverse=TRUE,ci.col='gray80',
                          line.col='black',lwd=1,plot=TRUE,
                          exterr=TRUE,model=1,show.ellipses=1*(model!=2),
                          growth=FALSE,hide=NULL,omit=NULL,omit.col=NA,...){
    y <- data2york(x,inverse=inverse)
    d2calc <- clear(y,hide,omit)
    fit <- regression(d2calc,model=model)
    out <- isochron_init(fit,alpha=alpha)
    out$y0[c('y','s[y]')] <- out$a
    if (inverse){
        R76 <- out$a
        x.lab <- quote(''^204*'Pb/'^206*'Pb')
        y.lab <- quote(''^207*'Pb/'^206*'Pb')
        out$y0label <- quote(''^207*'Pb/'^206*'Pb = ')
    } else {
        R76 <- out$b
        x.lab <- quote(''^206*'Pb/'^204*'Pb')
        y.lab <- quote(''^207*'Pb/'^204*'Pb')
        out$y0label <- quote('('^207*'Pb/'^204*'Pb)'[o]*' = ')
    }
    out$displabel <- quote('dispersion = ')
    out$age[c('t','s[t]')] <-
        get.Pb207Pb206.age(R76[1],R76[2],exterr=exterr)
    out <- ci_isochron(out)
    if (model==1){
        out$age['disp[t]'] <-
            out$fact*get.Pb207Pb206.age(R76[1],sqrt(out$mswd)*R76[2],
                                        exterr=exterr)[2]
    }
    if (plot) {
        scatterplot(y,xlim=xlim,ylim=ylim,alpha=alpha,
                    show.ellipses=show.ellipses,
                    show.numbers=show.numbers,levels=levels,
                    clabel=clabel,ellipse.col=ellipse.col,fit=out,
                    ci.col=ci.col,line.col=line.col,lwd=lwd,
                    hide=hide,omit=omit,omit.col=omit.col,...)
        if (growth) plot_PbPb_evolution(out$age['t'],inverse=inverse)
        graphics::title(isochrontitle(out,sigdig=sigdig,type='Pb-Pb'),
                        xlab=x.lab,ylab=y.lab)
    }
    invisible(out)
}
#' @rdname isochron
#' @export
isochron.RbSr <- function(x,xlim=NA,ylim=NA,alpha=0.05,sigdig=2,
                          show.numbers=FALSE,levels=NA,clabel="",
                          ellipse.col=c("#00FF0080","#FF000080"),
                          ci.col='gray80',line.col='black',lwd=1,
                          plot=TRUE,exterr=TRUE,model=1,
                          show.ellipses=1*(model!=2),hide=NULL,
                          omit=NULL,omit.col=NA,...){
    isochron_PD(x,'Rb87',xlim=xlim,ylim=ylim,alpha=alpha,
                sigdig=sigdig,show.numbers=show.numbers,
                levels=levels,clabel=clabel,ellipse.col=ellipse.col,
                ci.col=ci.col,line.col=line.col,lwd=lwd,plot=plot,
                exterr=exterr,model=model,show.ellipses=show.ellipses,
                hide=hide,omit=omit,omit.col=omit.col,...)
}
#' @rdname isochron
#' @export
isochron.ReOs <- function(x,xlim=NA,ylim=NA,alpha=0.05,sigdig=2,
                          show.numbers=FALSE,levels=NA,clabel="",
                          ellipse.col=c("#00FF0080","#FF000080"),
                          ci.col='gray80',line.col='black',lwd=1,
                          plot=TRUE,exterr=TRUE,model=1,
                          show.ellipses=1*(model!=2),hide=NULL,
                          omit=NULL,omit.col=NA,...){
    isochron_PD(x,'Re187',xlim=xlim,ylim=ylim,alpha=alpha,
                sigdig=sigdig,show.numbers=show.numbers,
                levels=levels,clabel=clabel,ellipse.col=ellipse.col,
                ci.col=ci.col,line.col=line.col,lwd=lwd,plot=plot,
                exterr=exterr,model=model,show.ellipses=show.ellipses,
                hide=hide,omit=omit,omit.col=omit.col,...)
}
#' @rdname isochron
#' @export
isochron.SmNd <- function(x,xlim=NA,ylim=NA,alpha=0.05,sigdig=2,
                          show.numbers=FALSE,levels=NA,clabel="",
                          ellipse.col=c("#00FF0080","#FF000080"),
                          ci.col='gray80',line.col='black',lwd=1,
                          plot=TRUE,exterr=TRUE,model=1,
                          show.ellipses=1*(model!=2),hide=NULL,
                          omit=NULL,omit.col=NA,...){
    isochron_PD(x,'Sm147',xlim=xlim,ylim=ylim,alpha=alpha,
                sigdig=sigdig,show.numbers=show.numbers,
                levels=levels,clabel=clabel,ellipse.col=ellipse.col,
                ci.col=ci.col,line.col=line.col,lwd=lwd,
                plot=plot,exterr=exterr,model=model,
                show.ellipses=show.ellipses,hide=hide,
                omit=omit,omit.col=omit.col,...)
}
#' @rdname isochron
#' @export
isochron.LuHf <- function(x,xlim=NA,ylim=NA,alpha=0.05,sigdig=2,
                          show.numbers=FALSE,levels=NA,clabel="",
                          ellipse.col=c("#00FF0080","#FF000080"),
                          ci.col='gray80',line.col='black',lwd=1,
                          plot=TRUE,exterr=TRUE,model=1,
                          show.ellipses=1*(model!=2),hide=NULL,
                          omit=NULL,omit.col=NA,...){
    isochron_PD(x,'Lu176',xlim=xlim,ylim=ylim,alpha=alpha,
                sigdig=sigdig,show.numbers=show.numbers,
                levels=levels,clabel=clabel,ellipse.col=ellipse.col,
                ci.col=ci.col,line.col=line.col,lwd=lwd,
                plot=plot,exterr=exterr,model=model,
                show.ellipses=show.ellipses,hide=hide,
                omit=omit,omit.col=omit.col,...)
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
                          clabel="",ellipse.col=c("#00FF0080","#FF000080"),
                          ci.col='gray80',line.col='black',lwd=1,
                          plot=TRUE,exterr=TRUE,model=1,
                          show.ellipses=1*(model!=2),hide=NULL,
                          omit=NULL,omit.col=NA,...){
    if (x$format %in% c(1,2)){
        out <- isochron_ThU_3D(x,type=type,model=model,
                               exterr=exterr,alpha=alpha,
                               hide=hide,omit=omit)
        intercept.type <- 'Th-U-3D'
    } else if (x$format %in% c(3,4)){
        out <- isochron_ThU_2D(x,type=type,model=model,
                               exterr=exterr,alpha=alpha,
                               hide=hide,omit=omit)
        intercept.type <- 'Th-U-2D'
    }
    if (type %in% c(1,3)){
        out$displabel <- quote('('^234*'U/'^232*'Th)-dispersion = ')
    } else if (type %in% c(2,4)){
        out$displabel <- quote('('^234*'U/'^238*'U)-dispersion = ')
    }
    if (plot){
        scatterplot(out$xyz,xlim=xlim,ylim=ylim,alpha=alpha,
                    show.numbers=show.numbers,levels=levels,
                    clabel=clabel,ellipse.col=ellipse.col,fit=out,
                    show.ellipses=show.ellipses,ci.col=ci.col,
                    line.col=line.col,lwd=lwd,
                    hide=hide,omit=omit,omit.col=omit.col,...)
        graphics::title(isochrontitle(out,sigdig=sigdig,
                        type=intercept.type,units='ka'),
                        xlab=out$xlab,ylab=out$ylab)
    }
    invisible(out)
}
#' @rdname isochron
#' @export
isochron.UThHe <- function(x,xlim=NA,ylim=NA,alpha=0.05,sigdig=2,
                           show.numbers=FALSE,levels=NA,clabel="",
                           ellipse.col=c("#00FF0080","#FF000080"),
                           ci.col='gray80',line.col='black',lwd=1,
                           plot=TRUE,model=1,show.ellipses=2*(model!=2),
                           hide=NULL,omit=NULL,omit.col='grey',...){
    y <- data2york(x)
    d2calc <- clear(y,hide,omit)
    fit <- regression(d2calc,model=model)
    out <- isochron_init(fit,alpha=alpha)
    out$y0[c('y','s[y]')] <- out$a
    out$age[c('t','s[t]')] <- out$b
    out <- ci_isochron(out)
    if (model==1)
        out$age['disp[t]'] <- out$fact*sqrt(out$mswd)*out$age['s[t]']
    out$displabel <- quote('He-dispersion = ')
    out$y0label <- quote('He'[o]*' = ')
    if (plot) {
        scatterplot(y,xlim=xlim,ylim=ylim,alpha=alpha,
                    show.numbers=show.numbers,levels=levels,
                    clabel=clabel,ellipse.col=ellipse.col,fit=out,
                    show.ellipses=show.ellipses,ci.col=ci.col,
                    line.col=line.col,lwd=lwd,hide=hide,
                    omit=omit,omit.col=omit.col,...)
        graphics::title(isochrontitle(out,sigdig=sigdig,type='U-Th-He'),
                        xlab="P",ylab="He")
    }
    invisible(out)
}

isochron_ThU_3D <- function(x,type=2,model=1,exterr=TRUE,
                            alpha=0.05,hide=NULL,omit=NULL){
    if (type == 1){ # 0/2 vs 8/2
        osmond <- FALSE
        ia <- 'A'
        ib <- 'B'
        i48 <- 'b'
        i08 <- 'B'
        id <- c('X','sX','Z','sZ','rXZ')
        x.lab <- quote(''^238*'U/'^232*'Th')
        y.lab <- quote(''^230*'Th/'^232*'Th')
    } else if (type == 2){ # 0/8 vs 2/8
        osmond <- TRUE
        ia <- 'A'
        ib <- 'B'
        i48 <- 'a'
        i08 <- 'A'
        id <- c('X','sX','Z','sZ','rXZ')
        x.lab <- quote(''^232*'Th/'^238*'U')
        y.lab <- quote(''^230*'Th/'^238*'U')
    } else if (type == 3){ # 4/2 vs 8/2
        osmond <- FALSE
        ia <- 'a'
        ib <- 'b'
        i48 <- 'b'
        i08 <- 'B'
        id <- c('X','sX','Y','sY','rXY')
        x.lab <- quote(''^238*'U/'^232*'Th')
        y.lab <- quote(''^234*'U/'^232*'Th')
    } else if (type == 4){ # 4/8 vs 2/8
        osmond <- TRUE
        ia <- 'a'
        ib <- 'b'
        i48 <- 'a'
        i08 <- 'A'
        id <- c('X','sX','Y','sY','rXY')
        x.lab <- quote(''^232*'Th/'^238*'U')
        y.lab <- quote(''^234*'U/'^238*'U')
    }
    tit <- data2tit(x,osmond=osmond)
    d2calc <- clear(tit,hide,omit)
    fit <- regression(d2calc,model=model,type="titterington")
    out <- isochron_init(fit,alpha=alpha)
    out$xyz <- tit
    out$a <- c(out$par[ia],sqrt(out$cov[ia,ia]))
    out$b <- c(out$par[ib],sqrt(out$cov[ib,ib]))
    out$cov.ab <- out$cov[ia,ib]
    tst <- get.ThU.age(out$par[i08],sqrt(out$cov[i08,i08]),
                       out$par[i48],sqrt(out$cov[i48,i48]),
                       out$cov[i48,i08],exterr=exterr)
    out$age['t'] <- tst['t']
    out$y0['y'] <- tst['48_0']
    out$age['s[t]'] <- tst['s[t]']
    out$y0['s[y]'] <- tst['s[48_0]']
    out$y0label <- quote('('^234*'U/'^238*'U)'[o]*'=')
    out <- ci_isochron(out,disp=FALSE)
    if (model==1 && out$mswd>1){
        tdispt <- get.ThU.age(out$par[i08],
                              sqrt(out$mswd)*sqrt(out$cov[i08,i08]),
                              out$par[i48],
                              sqrt(out$mswd)*sqrt(out$cov[i48,i48]),
                              out$mswd*out$cov[i48,i08],
                              exterr=exterr)
        out$age['disp[t]'] <- out$fact*tdispt['s[t]']
        out$y0['disp[y]'] <- out$fact*tdispt['s[48_0]']
    }
    out$xlab <- x.lab
    out$ylab <- y.lab
    out$xyz <- subset(out$xyz,select=id)
    out
}
isochron_ThU_2D <- function(x,type=2,model=1,exterr=TRUE,
                            alpha=0.05,hide=NULL,omit=NULL){
    y <- data2york(x,type=type)
    d2calc <- clear(y,hide,omit)
    fit <- regression(d2calc,model=model,type="york")
    out <- isochron_init(fit,alpha=alpha)
    out$xyz <- y
    if (type==1){
        Th230U238 <- out$b
        Th230Th232 <- out$a
        x.lab <- quote(''^238*'U/'^232*'Th')
        y.lab <- quote(''^230*'Th/'^232*'Th')
    } else if (type==2) {
        Th230U238 <- out$a
        Th230Th232 <- out$b
        x.lab <- quote(''^232*'Th/'^238*'U')
        y.lab <- quote(''^230*'Th/'^238*'U')
    }
    out$age[c('t','s[t]')] <-
        get.ThU.age(Th230U238[1],Th230U238[2],
                    exterr=exterr)[c('t','s[t]')]
    out$y0[c('y','s[y]')] <-
        get.Th230Th232_0x(out$age['t'],Th230Th232[1],Th230Th232[2])
    out <- ci_isochron(out,disp=FALSE)
    if (model==1 && out$mswd>1){
        out$age['disp[t]'] <-
            out$fact*get.ThU.age(Th230U238[1],
                                  sqrt(out$mswd)*Th230U238[2],
                                  exterr=exterr)['s[t]']
        out$y0['disp[y]'] <-
            out$fact*get.Th230Th232_0x(out$age['t'],Th230Th232[1],
                                        sqrt(out$mswd)*Th230Th232[2])[2]
    }
    out$xlab <- x.lab
    out$ylab <- y.lab
    out
}

isochron_PD <- function(x,nuclide,xlim=NA,ylim=NA,alpha=0.05,
                        sigdig=2,show.numbers=FALSE,levels=NA,
                        clabel="",ellipse.col=c("#00FF0080","#FF000080"),
                        ci.col='gray80',line.col='black',lwd=1,
                        plot=TRUE,exterr=TRUE,model=1,
                        show.ellipses=1*(model!=2),bratio=1,
                        hide=NULL,omit=NULL,...){
    y <- data2york(x)
    d2calc <- clear(y,hide,omit)
    fit <- regression(d2calc,model=model)
    out <- isochron_init(fit,alpha=alpha)
    out$y0[c('y','s[y]')] <- out$a
    out$age[c('t','s[t]')] <- get.PD.age(out$b['b'],out$b['s[b]'],
                                         nuclide,exterr=exterr,
                                         bratio=bratio)
    out <- ci_isochron(out)
    if (model==1){
        out$age['disp[t]'] <- out$fact*get.PD.age(out$b['b'],
                              sqrt(out$mswd)*out$b['s[b]'],
                              nuclide,exterr=exterr,bratio=bratio)[2]
    }
    if (identical(nuclide,'Sm147')){
        x.lab <- quote(''^147*'Sm/'^144*'Nd')
        y.lab <- quote(''^143*'Nd/'^144*'Nd')
    } else if (identical(nuclide,'Re187')){
        x.lab <- quote(''^187*'Re/'^188*'Os')
        y.lab <- quote(''^187*'Os/'^188*'Os')
    } else if (identical(nuclide,'Rb87')){
        x.lab <- quote(''^87*'Rb/'^86*'Sr')
        y.lab <- quote(''^87*'Sr/'^86*'Sr')
    } else if (identical(nuclide,'Lu176')){
        x.lab <- quote(''^176*'Lu/'^177*'Hf')
        y.lab <- quote(''^176*'Hf/'^177*'Hf')
    } else if (identical(nuclide,'K40')){
        x.lab <- quote(''^40*'K/'^44*'Ca')
        y.lab <- quote(''^40*'Ca/'^44*'Ca')
    }
    out$displabel <-
        substitute(a*b*c,list(a='(',b=y.lab,c=')-dispersion = '))
    out$y0label <-
        substitute(a*b*c,list(a='(',b=y.lab,c=quote(')'[o]*' = ')))
    if (plot){
        scatterplot(y,xlim=xlim,ylim=ylim,alpha=alpha,
                    show.ellipses=show.ellipses,
                    show.numbers=show.numbers,levels=levels,
                    clabel=clabel,ellipse.col=ellipse.col,fit=out,
                    ci.col=ci.col,line.col=line.col,lwd=lwd,
                    hide=hide,omit=omit,...)
        graphics::title(isochrontitle(out,sigdig=sigdig,type='PD'),
                        xlab=x.lab,ylab=y.lab)
    }
    invisible(out)
}

isochron_init <- function(fit,alpha=0.05){
    out <- fit
    if (fit$model==1){
        out$age <- rep(NA,4)
        out$y0 <- rep(NA,4)
        names(out$age) <- c('t','s[t]','ci[t]','disp[t]')
        names(out$y0) <- c('y','s[y]','ci[y]','disp[y]')
    } else {
        out$age <- rep(NA,3)
        out$y0 <- rep(NA,3)
        names(out$age) <- c('t','s[t]','ci[t]')
        names(out$y0) <- c('y','s[y]','ci[y]')
    }
    if (fit$model < 3){
        out$fact <- tfact(alpha,fit$df)
    } else {
        out$fact <- nfact(alpha)
        out$w <- c(out$w,NA,NA)
        names(out$w) <- c('s','ll','ul')
    }
    out$alpha <- alpha
    class(out) <- "isochron"
    out
}
regression_init <- function(fit,alpha=0.05){
    out <- fit
    out$displabel <- quote('y-dispersion = ')
    out$y0label <- quote('y-intercept = ')
    if (fit$model==1){
        out$a <- rep(NA,4)
        out$b <- rep(NA,4)
        names(out$a) <- c('a','s[a]','ci[a]','disp[a]')
        names(out$b) <- c('b','s[b]','ci[b]','disp[b]')
    } else {
        out$a <- rep(NA,3)
        out$b <- rep(NA,3)
        names(out$a) <- c('a','s[a]','ci[a]')
        names(out$b) <- c('b','s[b]','ci[b]')
    }
    if (fit$model < 3){
        out$fact <- tfact(alpha,fit$df)
    } else {
        out$fact <- nfact(alpha)
        out$w <- c(fit$w,NA,NA)
        names(out$w) <- c('s','ll','ul')
    }
    out$a[c('a','s[a]')] <- fit$a[c('a','s[a]')]
    out$b[c('b','s[b]')] <- fit$b[c('b','s[b]')]
    out$a['ci[a]'] <- out$fact*fit$a['s[a]']
    out$b['ci[b]'] <- out$fact*fit$b['s[b]']
    if (out$model==1){
        out$a['disp[a]'] <- out$fact*sqrt(fit$mswd)*fit$a['s[a]']
        out$b['disp[b]'] <- out$fact*sqrt(fit$mswd)*fit$b['s[b]']
    }
    out$alpha <- alpha
    class(out) <- "isochron"
    out
}

get.limits <- function(x,sx){
    minx <- min(x-3*sx,na.rm=TRUE)
    maxx <- max(x+3*sx,na.rm=TRUE)
    c(minx,maxx)
}

plot_PbPb_evolution <- function(tt,inverse=TRUE){
    nn <- 50
    tijd <- seq(0,tt,length.out=nn)
    ticks <- pretty(tijd)
    tijd[nn] <- max(ticks)
    xy <- stacey.kramers(tijd,inverse=inverse)
    graphics::lines(xy[,1],xy[,2])
    xy <- stacey.kramers(ticks,inverse=inverse)
    graphics::points(xy[,1],xy[,2],pch=20)
    graphics::text(xy[,1],xy[,2],labels=ticks,pos=3)
}

isochrontitle <- function(fit,sigdig=2,type=NA,units="Ma",...){
    if (fit$model==1 && fit$mswd>1){
        args1 <- quote(a%+-%b~'|'~c~'|'~d~u~'(n='*n*')')
        args2 <- quote(a%+-%b~'|'~c~'|'~d~u)
    } else {
        args1 <- quote(a%+-%b~'|'~c~u~'(n='*n*')')
        args2 <- quote(a%+-%b~'|'~c~u)
    }
    if (is.na(type)){
        intercept <- roundit(fit$a[1],fit$a[2:4],sigdig=sigdig)
        slope <- roundit(fit$b[1],fit$b[2:4],sigdig=sigdig)
        expr1 <- 'slope ='
        expr2 <- 'intercept ='
        list1 <- list(a=slope[1],
                      b=slope[2],
                      c=slope[3],
                      u='',
                      n=fit$n)
        list2 <- list(a=intercept[1],
                      b=intercept[2],
                      c=intercept[3],
                      u='')
        if (fit$model==1 && fit$mswd>1){
            list1$d <- slope[4]
            list2$d <- intercept[4]
        }
    } else {
        rounded.age <- roundit(fit$age[1],fit$age[2:4],sigdig=sigdig)
        rounded.intercept <- roundit(fit$y0[1],fit$y0[2:4],sigdig=sigdig)
        expr1 <- 'age ='
        list1 <- list(a=rounded.age[1],
                      b=rounded.age[2],
                      c=rounded.age[3],
                      u=units,
                      n=fit$n)
        list2 <- list(a=rounded.intercept[1],
                      b=rounded.intercept[2],
                      c=rounded.intercept[3],
                      u='')
        if (fit$model==1 && fit$mswd>1){
            list1$d <- rounded.age[4]
            list2$d <- rounded.intercept[4]
        }
        expr2 <- fit$y0label
    }
    call1 <- substitute(e~a,list(e=expr1,a=args1))
    call2 <- substitute(e~a,list(e=expr2,a=args2))
    line1 <- do.call(substitute,list(eval(call1),list1))
    line2 <- do.call(substitute,list(eval(call2),list2))
    if (fit$model==1){
        line3 <- substitute('MSWD ='~a~', p('*chi^2*')='~b,
                            list(a=signif(fit$mswd,sigdig),
                                 b=signif(fit$p.value,sigdig)))
        mymtext(line1,line=2,...)
        mymtext(line2,line=1,...)
        mymtext(line3,line=0,...)
    } else if (fit$model==2){
        mymtext(line1,line=1,...)
        mymtext(line2,line=0,...)
    } else if (fit$model==3){
        rounded.disp <- roundit(fit$w[1],fit$w[2:3],sigdig=sigdig)
        list3 <- list(a=rounded.disp[1],c=rounded.disp[2],b=rounded.disp[3])
        args3 <- quote(a+b/-c)
        expr3 <- fit$displabel
        call3 <- substitute(e~a,list(e=expr3,a=args3))
        line3 <- do.call(substitute,list(eval(call3),list3))
        mymtext(line1,line=2,...)
        mymtext(line2,line=1,...)
        mymtext(line3,line=0,...)
    }
}
