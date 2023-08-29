#' @title
#' Calculate and plot isochrons
#'
#' @description
#' Plots cogenetic U-Pb, Ar-Ar, K-Ca, Pb-Pb, Th-Pb, Rb-Sr, Sm-Nd,
#' Re-Os, Lu-Hf, U-Th-He or Th-U data as X-Y scatterplots, fits an
#' isochron curve through them using the \code{york},
#' \code{titterington} or \code{ludwig} function, and computes the
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
#' isochron.  The easiest of these is total least squares
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
#' Ar-Ar, K-Ca, Pb-Pb, Th-Pb, Rb-Sr, Sm-Nd, Re-Os and Lu-Hf
#' isochrons. The maximum likelihood algorithm of Titterington and
#' Halliday (1979) was generalised from two to three dimensions by
#' Ludwig and Titterington (1994) for U-series disequilibrium dating.
#' Also this algorithm is implemented in \code{IsoplotR}. Finally, the
#' constrained maximum likelihood algorithms of Ludwig (1998) and
#' Vermeesch (2020) are used for isochron regression of U-Pb data. The
#' extent to which the observed scatter in the data can be explained
#' by the analytical uncertainties can be assessed using the Mean
#' Square of the Weighted Deviates (MSWD, McIntyre et al., 1966),
#' which is defined as:
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
#' \item Ignore the analytical uncertainties and perform a total
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
#' \code{X}: the x-variable
#'
#' \code{sX}: the standard error of \code{X}
#'
#' \code{Y}: the y-variable
#'
#' \code{sY}: the standard error of \code{Y}
#'
#' \code{rXY}: the correlation coefficient of \code{X} and \code{Y}
#'
#' OR
#'
#' an object of class \code{ArAr}, \code{KCa}, \code{PbPb},
#' \code{UPb}, \code{ThPb}, \code{ReOs}, \code{RbSr}, \code{SmNd},
#' \code{LuHf}, \code{UThHe} or \code{ThU}.
#'
#' @param oerr indicates whether the analytical uncertainties of the
#'     output are reported in the plot title as:
#' 
#' \code{1}: 1\eqn{\sigma} absolute uncertainties.
#' 
#' \code{2}: 2\eqn{\sigma} absolute uncertainties.
#' 
#' \code{3}: absolute (1-\eqn{\alpha})\% confidence intervals, where
#' \eqn{\alpha} equales the value that is stored in
#' \code{settings('alpha')}.
#'
#' \code{4}: 1\eqn{\sigma} relative uncertainties (\eqn{\%}).
#' 
#' \code{5}: 2\eqn{\sigma} relative uncertainties (\eqn{\%}).
#'
#' \code{6}: relative (1-\eqn{\alpha})\% confidence intervals, where
#' \eqn{\alpha} equales the value that is stored in
#' \code{settings('alpha')}.
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
#' @param ellipse.fill
#' Fill colour for the error ellipses. This can either be a single
#' colour or multiple colours to form a colour ramp. Examples:
#'
#' a single colour: \code{rgb(0,1,0,0.5)}, \code{'#FF000080'},
#' \code{'white'}, etc.;
#'
#' multiple colours: \code{c(rbg(1,0,0,0.5)},
#' \code{rgb(0,1,0,0.5))}, \code{c('#FF000080','#00FF0080')},
#' \code{c('blue','red')}, \code{c('blue','yellow','red')}, etc.;
#'
#' a colour palette: \code{rainbow(n=100)},
#' \code{topo.colors(n=100,alpha=0.5)}, etc.; or
#'
#' a reversed palette: \code{rev(topo.colors(n=100,alpha=0.5))},
#' etc.
#'
#' For empty ellipses, set \code{ellipse.col=NA}
#' 
#' @param ellipse.stroke the stroke colour for the error
#'     ellipses. Follows the same formatting guidelines as
#'     \code{ellipse.fill}
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
#' \code{1}: Error-weighted least squares regression
#'
#' \code{2}: Total least squares regression
#'
#' \code{3}: Error-weighted least squares with overdispersion term
#'
#' @param wtype controls the parameter responsible for the
#'     overdispersion in model-3 regression.
#'
#' \code{0}, \code{'a'} or \code{'intercept'}: attributes the
#' overdispersion to the y-intercept of the isochron.
#'
#' \code{1}, \code{'b'} or \code{'slope'}: attributes the
#' overdispersion to the slope of the isochron.
#'
#' \code{'A'}: only available if \code{x} has class \code{ThU} and
#' \code{x$format} is 1 or 2. Attributes the overdispersion to the
#' authigenic \eqn{^{230}}Th/\eqn{^{238}}U-intercept of the isochron.
#'
#' \code{'B'}: only available if \code{x} has class \code{ThU} and
#' \code{x$format} is 1 or 2. Attributes the overdispersion to the
#' \eqn{^{230}}Th/\eqn{^{232}}Th-slope of the isochron.
#' 
#' @param show.ellipses show the data as:
#'
#' \code{0}: points
#'
#' \code{1}: error ellipses
#'
#' \code{2}: error crosses
#'
#' @param xlab text label for the horizontal plot axis
#' 
#' @param hide vector with indices of aliquots that should be removed
#'     from the plot.
#' 
#' @param omit vector with indices of aliquots that should be plotted
#'     but omitted from the isochron age calculation.
#' 
#' @param omit.fill fill colour that should be used for the omitted
#'     aliquots.
#' 
#' @param omit.stroke stroke colour that should be used for the
#'     omitted aliquots.
#'
#' @param ylab text label for the vertical plot axis
#' 
#' @param ... optional arguments to be passed on to the generic plot
#'     function if \code{model=2}
#'
#' @return If \code{x} has class \code{PbPb}, \code{ThPb},
#'     \code{ArAr}, \code{KCa}, \code{RbSr}, \code{SmNd}, \code{ReOs}
#'     or \code{LuHf}, or \code{UThHe}, returns a list with the
#'     following items:
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
#' \item{y0}{a two- or three-element list containing:
#'
#' \code{y}: the atmospheric \eqn{^{40}}Ar/\eqn{^{36}}Ar or initial
#' \eqn{^{40}}Ca/\eqn{^{44}}Ca, \eqn{^{187}}Os/\eqn{^{188}}Os,
#' \eqn{^{87}}Sr/\eqn{^{87}}Rb, \eqn{^{143}}Nd/\eqn{^{144}}Nd,
#' \eqn{^{176}}Hf/\eqn{^{177}}Hf or \eqn{^{208}}Pb/\eqn{^{204}}Pb
#' ratio.
#'
#' \code{s[y]}: the standard error of \code{y}
#'
#' \code{disp[y]}: the standard error of \code{y} enhanced by
#' \eqn{\sqrt{mswd}} (only applicable if \code{ model=1}).  }
#'
#' \item{age}{a three-element list containing:
#'
#' \code{t}: the \eqn{^{207}}Pb/\eqn{^{206}}Pb,
#' \eqn{^{208}}Pb/\eqn{^{232}}Th, \eqn{^{40}}Ar/\eqn{^{39}}Ar,
#' \eqn{^{40}}K/\eqn{^{40}}Ca, \eqn{^{187}}Os/\eqn{^{187}}Re,
#' \eqn{^{87}}Sr/\eqn{^{87}}Rb, \eqn{^{143}}Nd/\eqn{^{144}}Nd or
#' \eqn{^{176}}Hf/\eqn{^{177}}Hf age.
#'
#' \code{s[t]}: the standard error of \code{t}
#'
#' \code{disp[t]}: the standard error of \code{t} enhanced by
#' \eqn{\sqrt{mswd}} (only applicable if \code{ model=1}).  }
#'
#' \item{mswd}{the mean square of the residuals (a.k.a `reduced
#'     Chi-square') statistic (omitted if \code{model=2}).}
#'
#' \item{p.value}{the p-value of a Chi-square test for linearity
#' (omitted if \code{model=2})}
#'
#' \item{w}{the overdispersion term, i.e. a two-element vector with
#' the standard deviation of the (assumed) Normally distributed
#' geological scatter that underlies the measurements, and its
#' standard error (only returned if \code{model=3}).}
#'
#' \item{ski}{(only reported if \code{x} has class \code{PbPb} and
#' \code{growth} is \code{TRUE}) the intercept(s) of the isochron with
#' the Stacey-Kramers mantle evolution curve.}
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
#' \item{y0}{a three-element vector containing:
#'
#' \code{y}: the initial \eqn{^{234}}U/\eqn{^{238}}U-ratio
#'
#' \code{s[y]}: the standard error of \code{y}
#'
#' \code{disp[y]}: the standard error of \code{y} enhanced by
#' \eqn{\sqrt{mswd}}.}
#'
#' \item{age}{a two (or three) element vector containing:
#'
#' \code{t}: the initial \eqn{^{234}}U/\eqn{^{238}}U-ratio
#'
#' \code{s[t]}: the standard error of \code{t}
#'
#' \code{disp[t]}: the standard error of \code{t} enhanced by
#' \eqn{\sqrt{mswd}} (only reported if \code{model=1}).}
#'
#' \item{w}{the overdispersion term, i.e. a two-element vector with
#' the standard deviation of the (assumedly) Normally distributed
#' geological scatter that underlies the measurements, and its
#' standard error.}
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
#' OR if \code{x} has class \code{UPb}:
#'
#' \describe{
#'
#' \item{par}{if \code{model=1} or \code{2}, a three element vector
#' containing the isochron age and the common Pb isotope ratios. If
#' \code{model=3}, adds a fourth element with the overdispersion
#' parameter \eqn{w}.}
#'
#' \item{cov}{the covariance matrix of \code{par}}
#'
#' \item{logpar}{the logarithm of \code{par}}
#'
#' \item{logcov}{the logarithm of \code{cov}}
#'
#' \item{n}{the number of analyses in the dataset}
#'
#' \item{df}{the degrees of freedom for the linear fit, i.e. \eqn{2n-3}}
#'
#' \item{a}{the y-intercept and its standard error}
#'
#' \item{b}{the isochron slope and its standard error}
#'
#' \item{cov.ab}{the covariance between \code{a} and \code{b}.}
#'
#' \item{mswd}{the mean square of the residuals (a.k.a `reduced
#'     Chi-square') statistic.}
#'
#' \item{p.value}{the p-value of a Chi-square test for linearity.}
#'
#' \item{y0}{a two or three-element vector containing:
#'
#' \code{y}: the initial \eqn{^{206}}Pb/\eqn{^{204}}Pb-ratio (if
#' \code{type=1} and \code{x$format=4,5} or \code{6});
#' \eqn{^{207}}Pb/\eqn{^{204}}Pb-ratio (if \code{type=2} and
#' \code{x$format=4,5} or \code{6});
#' \eqn{^{208}}Pb/\eqn{^{206}}Pb-ratio (if \code{type=1} and
#' \code{x$format=7} or \code{8}); 
#' \eqn{^{208}}Pb/\eqn{^{207}}Pb-ratio (if \code{type=2} and
#' \code{x$format=7} or \code{8});
#' \eqn{^{206}}Pb/\eqn{^{208}}Pb-ratio (if \code{type=3} and
#' \code{x$format=7} or \code{8}); or
#' \eqn{^{207}}Pb/\eqn{^{208}}Pb-ratio (if \code{type=4} and
#' \code{x$format=7} or \code{8}).
#'
#' \code{s[y]}: the standard error of \code{y}
#'
#' \code{disp[y]}: the standard error of \code{y} enhanced by
#' \eqn{\sqrt{mswd}} (only returned if \code{model=1})}
#'
#' \item{y0label}{the y-axis label of the isochron plot}
#'
#' \item{age}{a two (or three) element vector containing:
#'
#' \code{t}: the isochron age
#'
#' \code{s[t]}: the standard error of \code{t}
#'
#' \code{disp[t]}: the standard error of \code{t} enhanced by
#' \eqn{\sqrt{mswd}} (only reported if \code{model=1}).}
#'
#' \item{xlab}{the x-label of the isochron plot}
#'
#' \item{ylab}{the y-label of the isochron plot}
#'
#' }
#'
#' @examples
#' attach(examples)
#' isochron(RbSr)
#'
#' fit <- isochron(ArAr,inverse=FALSE,plot=FALSE)
#'
#' dev.new()
#' isochron(ThU,type=4)
#'
#' @seealso
#' \code{\link{york}},
#' \code{\link{titterington}},
#' \code{\link{ludwig}}
#'
#' @references
#' Ludwig, K.R. and Titterington, D.M., 1994. Calculation of
#' \eqn{^{230}}Th/U isochrons, ages, and errors. Geochimica et
#' Cosmochimica Acta, 58(22), pp.5031-5042.
#'
#' Ludwig, K.R., 1998. On the treatment of concordant uranium-lead
#' ages. Geochimica et Cosmochimica Acta, 62(4), pp.665-676.
#'
#' Nicolaysen, L.O., 1961. Graphic interpretation of discordant age
#' measurements on metamorphic rocks. Annals of the New York Academy
#' of Sciences, 91(1), pp.198-206.
#'
#' Titterington, D.M. and Halliday, A.N., 1979. On the fitting of
#' parallel isochrons and the method of maximum likelihood. Chemical
#' Geology, 26(3), pp.183-195.
#'
#' Vermeesch, P., 2020. Unifying the U-Pb and Th-Pb methods: joint
#' isochron regression and common Pb correction, Geochronology, 2,
#' 119-131.
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
isochron.default <- function(x,oerr=3,sigdig=2,show.numbers=FALSE,
                             levels=NA,clabel="",xlab='x',ylab='y',
                             ellipse.fill=c("#00FF0080","#FF000080"),
                             ellipse.stroke='black',ci.col='gray80',
                             line.col='black',lwd=1,plot=TRUE,
                             title=TRUE,model=1,show.ellipses=1*(model!=2),
                             hide=NULL,omit=NULL,omit.fill=NA,
                             omit.stroke='grey',...){
    d2calc <- clear(x,hide,omit)
    fit <- regression(data2york(d2calc),model=model)
    genericisochronplot(x=x,fit=fit,oerr=oerr,sigdig=sigdig,
                        show.numbers=show.numbers,levels=levels,clabel=clabel,
                        xlab=xlab,ylab=ylab,ellipse.fill=ellipse.fill,
                        ellipse.stroke=ellipse.stroke,ci.col=ci.col,
                        line.col=line.col,lwd=lwd,plot=plot,title=title,
                        model=model,show.ellipses=1*(model!=2),
                        hide=hide,omit=omit,omit.fill=omit.fill,
                        omit.stroke=omit.stroke,...)
}
#' @rdname isochron
#' @export
isochron.other <- function(x,oerr=3,sigdig=2,show.numbers=FALSE,
                           levels=NA,clabel="",xlab='x',ylab='y',
                           ellipse.fill=c("#00FF0080","#FF000080"),
                           ellipse.stroke='black',ci.col='gray80',
                           line.col='black',lwd=1,plot=TRUE,
                           title=TRUE,model=1,show.ellipses=1*(model!=2),
                           hide=NULL,omit=NULL,omit.fill=NA,
                           omit.stroke='grey',...){
    d2calc <- clear(x,hide,omit)
    if (x$format%in%c(4,5)){
        yd <- data2york(d2calc$x,format=d2calc$format)
        fit <- regression(yd,model=model)
    } else if (x$format==6){
        fit <- regression(d2calc$x,model=model,type='ogls')
    } else {
        stop("Invalid data format for isochron regression.")
    }
    genericisochronplot(x=x,fit=fit,oerr=oerr,sigdig=sigdig,
                        show.numbers=show.numbers,levels=levels,clabel=clabel,
                        xlab=xlab,ylab=ylab,ellipse.fill=ellipse.fill,
                        ellipse.stroke=ellipse.stroke,ci.col=ci.col,
                        line.col=line.col,lwd=lwd,plot=plot,title=title,
                        model=model,show.ellipses=1*(model!=2),
                        hide=hide,omit=omit,omit.fill=omit.fill,
                        omit.stroke=omit.stroke,...)
}
genericisochronplot <- function(x,fit,oerr=3,sigdig=2,show.numbers=FALSE,
                                levels=NA,clabel="",xlab='x',ylab='y',
                                ellipse.fill=c("#00FF0080","#FF000080"),
                                ellipse.stroke='black',ci.col='gray80',
                                line.col='black',lwd=1,plot=TRUE,
                                title=TRUE,show.ellipses=TRUE,
                                hide=NULL,omit=NULL,omit.fill=NA,
                                omit.stroke='grey',...){
    if (inflate(fit)){
        fit$a['disp[a]'] <- sqrt(fit$mswd)*fit$a['s[a]']
        fit$b['disp[b]'] <- sqrt(fit$mswd)*fit$b['s[b]']
    }
    if (plot){
        y <- data2york(x)
        scatterplot(y,oerr=oerr,show.ellipses=show.ellipses,
                    show.numbers=show.numbers,levels=levels,
                    clabel=clabel,ellipse.fill=ellipse.fill,
                    ellipse.stroke=ellipse.stroke,fit=fit,
                    ci.col=ci.col,line.col=line.col,lwd=lwd,
                    hide=hide,omit=omit,omit.fill=omit.fill,
                    omit.stroke=omit.stroke,...)
        if (title)
            graphics::title(isochrontitle(fit,oerr=oerr,sigdig=sigdig,units=''),
                            xlab=xlab,ylab=ylab)
    }
    invisible(fit)
}
#' @param anchor control parameters to fix the intercept age or common
#'     Pb composition of the isochron fit. This can be a scalar or a
#'     vector.
#'
#' If \code{anchor[1]=0}: do not anchor the isochron.
#'
#' If \code{anchor[1]=1}: fix the common Pb composition at the values
#' stored in \code{settings('iratio',...)}.
#'
#' If \code{anchor[1]=2}: force the isochron line to intersect the
#' concordia line at an age equal to \code{anchor[2]}.
#' 
#' @param type if \code{x} has class \code{UPb} and \code{x$format=4},
#'     \code{5} or \code{6}:
#'
#' \code{1}: \eqn{^{204}}Pb/\eqn{^{206}}Pb vs. \eqn{^{238}}U/\eqn{^{206}}Pb
#'
#' \code{2}: \eqn{^{204}}Pb/\eqn{^{207}}Pb vs. \eqn{^{235}}U/\eqn{^{207}}Pb
#'
#' if \code{x} has class \code{UPb} and \code{x$format=7} or \code{8}:
#'
#' \code{1}: \eqn{^{208}}Pb\eqn{{}_\circ}/\eqn{^{206}}Pb vs. \eqn{^{238}}U/\eqn{^{206}}Pb
#'
#' \code{2}: \eqn{^{208}}Pb\eqn{{}_\circ}/\eqn{^{207}}Pb vs. \eqn{^{235}}U/\eqn{^{207}}Pb
#' 
#' \code{3}: \eqn{^{206}}Pb\eqn{{}_\circ}/\eqn{^{208}}Pb
#' vs. \eqn{^{232}}Th/\eqn{^{208}}Pb
#'
#' \code{4}: \eqn{^{207}}Pb\eqn{{}_\circ}/\eqn{^{208}}Pb
#' vs. \eqn{^{232}}Th/\eqn{^{208}}Pb
#' 
#' if \code{x} has class \code{ThU}, and following the classification
#' of Ludwig and Titterington (1994), one of either:
#'
#' \code{1}: `Rosholt type-II' isochron, setting out
#' \eqn{^{230}}Th/\eqn{^{232}}Th vs. \eqn{^{238}}U/\eqn{^{232}}Th
#'
#' \code{2}: `Osmond type-II' isochron, setting out \eqn{^{230}}Th/\eqn{^{238}}U
#' vs. \eqn{^{232}}Th/\eqn{^{238}}U
#'
#' \code{3}: `Rosholt type-II' isochron, setting out \eqn{^{234}}U/\eqn{^{232}}Th
#' vs. \eqn{^{238}}U/\eqn{^{232}}Th
#'
#' \code{4}: `Osmond type-II' isochron, setting out
#' \eqn{^{234}}U/\eqn{^{238}}U vs. \eqn{^{232}}Th/\eqn{^{238}}U
#' 
#' @param joint logical. Only applies to U-Pb data formats 4 and
#'     above. If \code{TRUE}, carries out three dimensional
#'     regression.  If \code{FALSE}, uses two dimensional isochron
#'     regression.  The latter can be used to compute
#'     \eqn{{}^{207}}Pb/\eqn{{}^{235}}U isochrons, which are immune to
#'     the complexities of initial \eqn{{}^{234}}U/\eqn{{}^{238}}U
#'     disequilibrium.
#'
#' @param y0option controls the type of y-intercept or activity ratio
#'     that is reported along with the isochron age. Only relevant to
#'     U-Pb data and Th-U data formats 1 and 2.
#'
#' For U-Pb data:
#'
#' \code{y0option=1} reports the common Pb composition,
#'
#' \code{y0option=2} reports the initial \eqn{^{234}}U/\eqn{^{238}}U
#' activity ratio.
#'
#' \code{y0option=3} reports the initial \eqn{^{230}}Th/\eqn{^{238}}U
#' activity ratio,
#'
#' For Th-U data:
#'
#' \code{y0option=1} reports the authigenic
#' \eqn{^{234}}U/\eqn{^{238}}U activity ratio,
#'
#' \code{y0option=2} reports the detrital
#' \eqn{^{230}}Th/\eqn{^{232}}Th activity ratio,
#'
#' \code{y0option=3} reports the authigenic
#' \eqn{^{230}}Th/\eqn{^{238}}U activity ratio,
#' 
#' \code{y0option=4} reports the initial \eqn{^{234}}U/\eqn{^{238}}U
#' activity ratio.
#'
#' @rdname isochron
#' @export
isochron.UPb <- function(x,oerr=3,sigdig=2,show.numbers=FALSE,
                         levels=NA,clabel="",joint=TRUE,
                         ellipse.fill=c("#00FF0080","#FF000080"),
                         ellipse.stroke='black',type=1,
                         ci.col='gray80',line.col='black',lwd=1,
                         plot=TRUE,exterr=FALSE,model=1,
                         show.ellipses=1*(model!=2),anchor=0,
                         hide=NULL,omit=NULL,omit.fill=NA,
                         omit.stroke='grey',y0option=1,...){
    if (x$format<4){
        if (plot){
            out <- concordia_helper(x,type=2,show.age=model+1,oerr=oerr,
                                    sigdig=sigdig,show.numbers=show.numbers,
                                    levels=levels,clabel=clabel,
                                    ellipse.fill=ellipse.fill,
                                    ellipse.stroke=ellipse.stroke,exterr=exterr,
                                    anchor=anchor,hide=hide,omit=omit,
                                    y0option=y0option,omit.fill=omit.fill,
                                    omit.stroke=omit.stroke,...)
        } else {
            out <- ludwig(x,exterr=exterr,model=model,anchor=anchor)
        }
    } else {
        ns <- length(x)
        calcit <- (1:ns)%ni%c(hide,omit)
        x2calc <- subset(x,subset=calcit)
        fit <- ludwig(x2calc,model=model,anchor=anchor,
                      exterr=exterr,type=ifelse(joint,0,type))
        tt <- fit$par['t']
        a0 <- fit$par['a0']
        b0 <- fit$par['b0']
        l8 <- settings('lambda','U238')[1]
        l5 <- settings('lambda','U235')[1]
        md <- mediand(x$d)
        if (md$U48$option==2) md$U48 <- list(x=unname(fit$par['U48i']),option=1)
        if (md$ThU$option==2) md$ThU <- list(x=unname(fit$par['ThUi']),option=1)
        McL <- mclean(tt,d=md,exterr=exterr)
        if (type==1){                           # 04-08c/06 vs. 38/06
            x0inv <- McL$Pb206U238
            dx0invdt <- McL$dPb206U238dt
            E <- fit$cov[c('t','a0'),c('t','a0')]
            x.lab <- quote(''^238*'U/'^206*'Pb')
        } else if (type==2){                    # 04-08c/07 vs. 35/07
            x0inv <- McL$Pb207U235
            dx0invdt <- McL$dPb207U235dt
            E <- fit$cov[c('t','b0'),c('t','b0')]
            x.lab <- quote(''^235*'U/'^207*'Pb')
        } else if (type==3 & x$format%in%c(7,8)){  # 06c/08 vs. 32/08
            x0inv <- age_to_Pb208Th232_ratio(tt=tt,st=0)[1]
            dx0invdt <- McL$dPb208Th232dt
            E <- fit$cov[c('t','a0'),c('t','a0')]
            x.lab <- quote(''^232*'Th/'^208*'Pb')
        } else if (type==4 & x$format%in%c(7,8)){  # 07c/08 vs. 32/08
            x0inv <- age_to_Pb208Th232_ratio(tt=tt,st=0)[1]
            dx0invdt <- McL$dPb208Th232dt
            E <- fit$cov[c('t','b0'),c('t','b0')]
            x.lab <- quote(''^232*'Th/'^208*'Pb')
        } else {
            stop('Invalid isochron type.')
        }
        out <- fit
        out$age <- NULL
        out$age['t'] <- tt
        out$age['s[t]'] <- sqrt(fit$cov['t','t'])
        J <- matrix(0,2,2)
        if (x$format%in%c(4,5,6) & type==1){        # 04/06 vs. 38/06
            XY <- data2york(x,option=3)
            a <- 1/fit$par['a0']
            J[1,2] <- -a^2
            y.lab <- quote(''^204*'Pb/'^206*'Pb')
        } else if (x$format%in%c(4,5,6) & type==2){ # 04/07 vs. 35/07
            XY <- data2york(x,option=4)
            a <- 1/fit$par['b0']
            J[1,2] <- -a^2
            y.lab <- quote(''^204*'Pb/'^207*'Pb')
        } else if (x$format%in%c(7,8) & type==1){   # 08/06 vs. 38/06
            XY <- data2york(x,option=6,tt=tt)
            a <- 1/fit$par['a0']
            J[1,2] <- -a^2
            y.lab <- quote(''^208*'Pb'[c]*'/'^206*'Pb')
        } else if (x$format%in%c(7,8) & type==2){   # 08/07 vs. 35/07
            XY <- data2york(x,option=7,tt=tt)
            U <- settings('iratio','U238U235')[1]
            a <- 1/fit$par['b0']
            J[1,2] <- -a^2
            y.lab <- quote(''^208*'Pb'[c]*'/'^207*'Pb')
        } else if (x$format%in%c(7,8) & type==3){   # 06c/08 vs. 32/08
            XY <- data2york(x,option=8,tt=tt)
            a <- fit$par['a0']
            J[1,2] <- 1
            y.lab <- quote(''^206*'Pb'[c]*'/'^208*'Pb')
        } else if (x$format%in%c(7,8) & type==4){   # 07c/08 vs. 32/08
            XY <- data2york(x,option=9,tt=tt)
            a <- fit$par['b0']
            J[1,2] <- 1
            y.lab <- quote(''^207*'Pb'[c]*'/'^208*'Pb')
        } else {
            stop('Isochron regression is not available for this input format.')
        }
        out <- getUPby0(out=out,fmt=x$format,type=type,option=y0option)
        b <- -a*x0inv
        J[2,1] <- -a*dx0invdt
        J[2,2] <- x0inv*a^2
        cov.ab <- J%*%E%*%t(J)
        out$a <- c(a,sqrt(cov.ab[1,1]))
        out$b <- c(b,sqrt(cov.ab[2,2]))
        names(out$a) <- c('a','s[a]')
        names(out$b) <- c('b','s[b]')
        out$cov.ab <- cov.ab[1,2]
        if (inflate(out)){
            out$age['disp[t]'] <- sqrt(out$mswd)*out$age['s[t]']
            out$y0['disp[y]'] <- sqrt(out$mswd)*out$y0['s[y]']
        }
        if (plot){
            scatterplot(XY,oerr=oerr,show.ellipses=show.ellipses,
                        show.numbers=show.numbers,levels=levels,
                        clabel=clabel,ellipse.fill=ellipse.fill,
                        ellipse.stroke=ellipse.stroke,fit=out,
                        ci.col=ci.col,line.col=line.col,lwd=lwd,
                        hide=hide,omit=omit,omit.fill=omit.fill,
                        omit.stroke=omit.stroke,...)
            graphics::title(isochrontitle(out,oerr=oerr,sigdig=sigdig,type='U-Pb',
                                          y0option=y0option,dispunits=' Ma'),
                            xlab=x.lab,ylab=y.lab)
        }
    }
    invisible(out)
}

getUPby0 <- function(out,fmt=1,type=1,option=1){
    out$y0 <- c('y'=NA,'s[y]'=NA)
    if (option==1){
        if (fmt<4){                              # 07/06 vs. 38/06
            out$y0['y'] <- out$par['a0']
            out$y0['s[y]'] <- out$err['s','a0']
            out$y0label <- quote('('^207*'Pb/'^206*'Pb)'[c]*'=')
        } else if (fmt %in% c(4,5,6) & type==1){ # 04/06 vs. 38/06
            out$y0['y'] <- out$par['a0']
            out$y0['s[y]'] <- sqrt(out$cov['a0','a0'])
            out$y0label <- quote('('^206*'Pb/'^204*'Pb)'[c]*'=')
        } else if (fmt %in% c(4,5,6) & type==2){ # 04/07 vs. 35/07
            out$y0['y'] <- out$par['b0']
            out$y0['s[y]'] <- sqrt(out$cov['b0','b0'])
            out$y0label <- quote('('^207*'Pb/'^204*'Pb)'[c]*'=')
        } else if (fmt %in% c(7,8) & type==1){   # 08/06 vs. 38/06
            out$y0['y'] <- 1/out$par['a0']
            out$y0['s[y]'] <- out$y0[1]*sqrt(out$cov['a0','a0'])/out$par['a0']
            out$y0label <- quote('('^208*'Pb/'^206*'Pb)'[c]*'=')
        } else if (fmt %in% c(7,8) & type==2){   # 08/07 vs. 35/07
            out$y0['y'] <- 1/out$par['b0']
            out$y0['s[y]'] <- out$y0[1]*sqrt(out$cov['b0','b0'])/out$par['b0']
            out$y0label <- quote('('^208*'Pb/'^207*'Pb)'[c]*'=')
        } else if (fmt %in% c(7,8) & type==3){   # 06c/08 vs. 32/08
            out$y0['y'] <- out$par['a0']
            out$y0['s[y]'] <- sqrt(out$cov['a0','a0'])
            out$y0label <- quote('('^206*'Pb/'^208*'Pb)'[c]*'=')
        } else if (fmt %in% c(7,8) & type==4){   # 07c/08 vs. 32/08
            out$y0['y'] <- out$par['b0']
            out$y0['s[y]'] <- sqrt(out$cov['b0','b0'])
            out$y0label <- quote('('^207*'Pb/'^208*'Pb)'[c]*'=')
        }
    } else {
        if (option==2){
            out$y0label <- quote('('^234*'U/'^238*'U)'[i]*'=')
            y0par <- 'U48i'
        } else if (option==3){
            out$y0label <- quote('('^230*'Th/'^238*'U)'[i]*'=')
            y0par <- 'ThUi'
        } else {
            stop('Invalid y0option value.')
        }
        if (y0par %in% names(out$par)){
            out$y0['y'] <- out$par[y0par][1]
            out$y0['s[y]'] <- sqrt(out$cov[y0par,y0par])
        } else {
            out$y0['y'] <- 1
            out$y0['s[y]'] <- 0
        }
    }
    if (inflate(out)){
        out$y0['disp'] <- sqrt(out$mswd)*out$y0['s[y]']
    }
    out
}
#' @param inverse toggles between normal and inverse isochrons. If the
#'     isochron plots \code{Y} against \code{X}, and
#'
#' If \code{inverse=TRUE}, then \code{X} =
#' \eqn{{}^{204}}Pb/\eqn{{}^{206}}Pb and \code{Y} =
#' \eqn{{}^{207}}Pb/\eqn{{}^{206}}Pb (if \code{x} has class
#' \code{PbPb}), or \code{X} = \eqn{{}^{232}}Th/\eqn{{}^{208}}Pb and
#' \code{Y} = \eqn{{}^{204}}Pb/\eqn{{}^{208}}Pb (if \code{x} has class
#' \code{ThPb}), or \code{X} = \eqn{{}^{39}}Ar/\eqn{{}^{40}}Ar and
#' \code{Y} = \eqn{{}^{36}}Ar/\eqn{{}^{40}}Ar (if \code{x} has class
#' \code{ArAr}), or \code{X} = \eqn{{}^{40}}K/\eqn{{}^{40}}Ca and
#' \code{Y} = \eqn{{}^{44}}Ca/\eqn{{}^{40}}Ca (if \code{x} has class
#' \code{KCa}), or \code{X} = \eqn{{}^{87}}Rb/\eqn{{}^{87}}Sr and
#' \code{Y} = \eqn{{}^{86}}Sr/\eqn{{}^{87}}Sr (if \code{x} has class
#' \code{RbSr}), or \code{X} = \eqn{{}^{147}}Sm/\eqn{{}^{143}}Nd and
#' \code{Y} = \eqn{{}^{144}}Nd/\eqn{{}^{143}}Nd (if \code{x} has class
#' \code{SmNd}), or \code{X} = \eqn{{}^{187}}Re/\eqn{{}^{187}}Os and
#' \code{Y} = \eqn{{}^{188}}Os/\eqn{{}^{187}}Os (if \code{x} has class
#' \code{ReOs}), or \code{X} = \eqn{{}^{176}}Lu/\eqn{{}^{176}}Hf and
#' \code{Y} = \eqn{{}^{177}}Hf/\eqn{{}^{176}}Hf (if \code{x} has class
#' \code{LuHf}).
#' 
#' If \code{inverse=FALSE}, then \code{X} =
#' \eqn{{}^{206}}Pb/\eqn{{}^{204}}Pb and \code{Y} =
#' \eqn{{}^{207}}Pb/\eqn{{}^{204}}Pb (if \code{x} has class
#' \code{PbPb}), or \code{X} = \eqn{{}^{232}}Th/\eqn{{}^{204}}Pb and
#' \code{Y} = \eqn{{}^{208}}Pb/\eqn{{}^{204}}Pb (if \code{x} has class
#' \code{ThPb}), or \code{X} = \eqn{{}^{39}}Ar/\eqn{{}^{36}}Ar and
#' \code{Y} = \eqn{{}^{40}}Ar/\eqn{{}^{36}}Ar (if \code{x} has class
#' \code{ArAr}), or \code{X} = \eqn{{}^{40}}K/\eqn{{}^{44}}Ca and
#' \code{Y} = \eqn{{}^{40}}Ca/\eqn{{}^{44}}Ca (if \code{x} has class
#' \code{KCa}), or \code{X} = \eqn{{}^{87}}Rb/\eqn{{}^{86}}Sr and
#' \code{Y} = \eqn{{}^{87}}Sr/\eqn{{}^{86}}Sr (if \code{x} has class
#' \code{RbSr}), or \code{X} = \eqn{{}^{147}}Sm/\eqn{{}^{144}}Nd and
#' \code{Y} = \eqn{{}^{143}}Nd/\eqn{{}^{144}}Nd (if \code{x} has class
#' \code{SmNd}), or \code{X} = \eqn{{}^{187}}Re/\eqn{{}^{188}}Os and
#' \code{Y} = \eqn{{}^{187}}Os/\eqn{{}^{188}}Os (if \code{x} has class
#' \code{ReOs}), or \code{X} = \eqn{{}^{176}}Lu/\eqn{{}^{177}}Hf and
#' \code{Y} = \eqn{{}^{176}}Hf/\eqn{{}^{177}}Hf (if \code{x} has class
#' \code{LuHf}).
#'
#' @param exterr propagate external sources of uncertainty
#' (J, decay constant)?
#' 
#' @param growth add Stacey-Kramers Pb-evolution curve to the plot?
#' @rdname isochron
#' @export
isochron.PbPb <- function(x,oerr=3,sigdig=2,show.numbers=FALSE,levels=NA,
                          clabel="",ellipse.fill=c("#00FF0080","#FF000080"),
                          ellipse.stroke='black',inverse=TRUE,
                          ci.col='gray80',line.col='black',lwd=1,
                          plot=TRUE,exterr=TRUE,model=1,wtype='intercept',
                          growth=FALSE,show.ellipses=1*(model!=2),hide=NULL,
                          omit=NULL,omit.fill=NA,omit.stroke='grey',...){
    y <- data2york(x,inverse=inverse)
    d2calc <- clear(y,hide,omit)
    out <- regression(d2calc,model=model,wtype=wtype)
    if (inverse){
        R76 <- out$a
        out$y0[c('y','s[y]')] <- out$b
        x.lab <- quote(''^204*'Pb/'^206*'Pb')
        y.lab <- quote(''^207*'Pb/'^206*'Pb')
    } else {
        R76 <- out$b
        out$y0[c('y','s[y]')] <- out$a
        x.lab <- quote(''^206*'Pb/'^204*'Pb')
        y.lab <- quote(''^207*'Pb/'^204*'Pb')
    }
    out$y0label <- quote('('^207*'Pb/'^204*'Pb)'[c]*'=')
    out$age[c('t','s[t]')] <- get.Pb207Pb206.age(R76[1],R76[2],exterr=exterr)
    if (inflate(out)){
        out$age['disp[t]'] <- 
            get.Pb207Pb206.age(R76[1],sqrt(out$mswd)*R76[2],exterr=exterr)[2]
        out$y0['disp[y]'] <- sqrt(out$mswd)*out$y0['s[y]']
    }
    disp2age <- (inverse && wtype%in%c('intercept',0,'a')) ||
                (!inverse && wtype%in%c('slope',1,'b'))
    if (model==3 && disp2age){
        out$disp <- out$disp/mclean(out$age['t'])$dPb207Pb206dt
        dispunits <- ' Ma'
    } else {
        dispunits <- ''   
    }
    if (plot) {
        scatterplot(y,oerr=oerr,show.ellipses=show.ellipses,
                    show.numbers=show.numbers,levels=levels,
                    clabel=clabel,ellipse.fill=ellipse.fill,
                    ellipse.stroke=ellipse.stroke,fit=out,
                    ci.col=ci.col,line.col=line.col,lwd=lwd,
                    hide=hide,omit=omit,omit.fill=omit.fill,
                    omit.stroke=omit.stroke,...)
        if (growth){
            xylim <- graphics::par('usr')
            if (xylim[1]<0) xylim[1] <- xylim[2]/100
            if (xylim[3]<0) xylim[3] <- xylim[4]/100
            if (inverse){
                Pb64 <- 1/xylim[1:2]
                Pb74 <- xylim[3:4]/xylim[2:1]
            } else {
                Pb64 <- xylim[1:2]
                Pb74 <- xylim[3:4]
            }
            tx <- sk2t(Pb206Pb204=Pb64)
            ty <- sk2t(Pb207Pb204=Pb74)
            tmin <- max(min(tx),min(ty))
            tmax <- min(max(tx),max(ty))
            plot_PbPb_evolution(from=tmin,to=tmax,inverse=inverse)
            out$ski <- SK.intersection(out,inverse=inverse)
        } else {
            out$ski <- NULL
        }
        tit <- isochrontitle(out,oerr=oerr,sigdig=sigdig,type='Pb-Pb',
                             dispunits=dispunits,ski=out$ski)
        graphics::title(tit,xlab=x.lab,ylab=y.lab)
    }
    invisible(out)
}
SK.intersection <- function(fit,inverse,m=0,M=5000){
    SKi.misfit <- function(tt,a,b){
        i6474 <- stacey.kramers(tt)
        pred74 <- a + b*i6474[1]
        pred74-i6474[2]
    }
    if (inverse){
        a <- fit$b[1]
        b <- fit$y0[1]
    } else {
        a <- fit$a[1]
        b <- fit$b[1]
    }
    if ((M-m)<10){
        out <- NULL
    } else if (sign(SKi.misfit(m,a,b)) == sign(SKi.misfit(M,a,b))){
        ski1 <- SK.intersection(fit,inverse,m=m,M=m+(M-m)/2)
        ski2 <- SK.intersection(fit,inverse,m=m+(M-m)/2,M=M)
        out <- unique(c(ski1,ski2))
    } else {
        out <- stats::uniroot(SKi.misfit,interval=c(m,M),a=a,b=b)$root
    }
    out
}
#' @rdname isochron
#' @export
isochron.ArAr <- function(x,oerr=3,sigdig=2,show.numbers=FALSE,levels=NA,
                          clabel="",ellipse.fill=c("#00FF0080","#FF000080"),
                          ellipse.stroke='black',inverse=TRUE,
                          ci.col='gray80',line.col='black',lwd=1,
                          plot=TRUE,exterr=TRUE,model=1,wtype='intercept',
                          show.ellipses=1*(model!=2),hide=NULL,
                          omit=NULL,omit.fill=NA,omit.stroke='grey',...){
    y <- data2york(x,inverse=inverse)
    d2calc <- clear(y,hide,omit)
    out <- regression(d2calc,model=model,wtype=wtype)
    a <- out$a['a']
    sa <- out$a['s[a]']
    b <- out$b['b']
    sb <- out$b['s[b]']
    if (inverse) {
        R09 <- -b/a
        sR09 <- R09*sqrt((sa/a)^2 + (sb/b)^2 - 2*out$cov.ab/(a*b))
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
    out$y0label <- quote('('^40*'Ar/'^36*'Ar)'[0]*'=')
    out$age[c('t','s[t]')] <- get.ArAr.age(R09,sR09,x$J[1],x$J[2],exterr=exterr)
    if (inflate(out)){
        out$age['disp[t]'] <- get.ArAr.age(R09,sqrt(out$mswd)*sR09,
                                           x$J[1],x$J[2],exterr=exterr)[2]
        out$y0['disp[y]'] <- sqrt(out$mswd)*out$y0['s[y]']
    }
    dispunits <- ''
    if (model==3){
        if (wtype%in%c('slope',1,'b')){
            l40 <- lambda('K40')[1]
            dtd09 <- (x$J[1]/l40)/(x$J[1]*R09+1)
            d09db <- ifelse(inverse,1/a,1)
            out$disp <- dtd09*d09db*out$disp
            dispunits <- ' Ma'
        } else if (inverse){ # wtype%in%c('intercept',0,'a')
            w <- out$disp[1]
            sw <- out$disp[2]
            out$disp[1] <- w/a^2
            out$disp[2] <- out$disp[1]*sw/w
        }
    }
    if (plot) {
        scatterplot(y,oerr=oerr,show.ellipses=show.ellipses,
                    show.numbers=show.numbers,levels=levels,
                    clabel=clabel,ellipse.fill=ellipse.fill,
                    ellipse.stroke=ellipse.stroke,fit=out,
                    ci.col=ci.col,line.col=line.col,lwd=lwd,
                    hide=hide,omit=omit,omit.fill=omit.fill,
                    omit.stroke=omit.stroke,...)
        graphics::title(isochrontitle(out,oerr=oerr,sigdig=sigdig,
                                      dispunits=dispunits,type='Ar-Ar'),
                        xlab=x.lab,ylab=y.lab)
    }
    invisible(out)
}
#' @rdname isochron
#' @export
isochron.ThPb <- function(x,oerr=3,sigdig=2,show.numbers=FALSE,levels=NA,
                          clabel="",ellipse.fill=c("#00FF0080","#FF000080"),
                          ellipse.stroke='black',inverse=FALSE,
                          ci.col='gray80',line.col='black',lwd=1,
                          plot=TRUE,exterr=TRUE,model=1,wtype='intercept',
                          show.ellipses=1*(model!=2),hide=NULL,
                          omit=NULL,omit.fill=NA,omit.stroke='grey',...){
    isochron_PD(x,nuclide='Th232',oerr=oerr,sigdig=sigdig,
                show.numbers=show.numbers,levels=levels,
                clabel=clabel,ellipse.fill=ellipse.fill,
                ellipse.stroke=ellipse.stroke,inverse=inverse,
                ci.col=ci.col,line.col=line.col,lwd=lwd,plot=plot,
                exterr=exterr,model=model,wtype=wtype,
                show.ellipses=show.ellipses,hide=hide,omit=omit,
                omit.fill=omit.fill,omit.stroke=omit.stroke,...)
}
#' @rdname isochron
#' @export
isochron.KCa <- function(x,oerr=3,sigdig=2,show.numbers=FALSE,levels=NA,
                         clabel="",inverse=FALSE,ci.col='gray80',
                         ellipse.fill=c("#00FF0080","#FF000080"),
                         ellipse.stroke='black',line.col='black',
                         lwd=1, plot=TRUE,exterr=TRUE,model=1,wtype='intercept',
                         show.ellipses=1*(model!=2),hide=NULL,
                         omit=NULL,omit.fill=NA,omit.stroke='grey',...){
    isochron_PD(x,nuclide='K40',oerr=oerr,sigdig=sigdig,
                show.numbers=show.numbers,levels=levels,
                clabel=clabel,ellipse.fill=ellipse.fill,
                ellipse.stroke=ellipse.stroke,inverse=inverse,
                ci.col=ci.col,line.col=line.col,lwd=lwd,plot=plot,
                exterr=exterr,model=model,wtype=wtype,
                show.ellipses=show.ellipses,bratio=0.895,hide=hide,
                omit=omit,omit.fill=omit.fill,omit.stroke=omit.stroke,...)
}
#' @rdname isochron
#' @export
isochron.RbSr <- function(x,oerr=3,sigdig=2,show.numbers=FALSE,levels=NA,
                          clabel="",ellipse.fill=c("#00FF0080","#FF000080"),
                          ellipse.stroke='black',inverse=FALSE,
                          ci.col='gray80',line.col='black',
                          lwd=1,plot=TRUE,exterr=TRUE,model=1,wtype='intercept',
                          show.ellipses=1*(model!=2),hide=NULL,
                          omit=NULL,omit.fill=NA,omit.stroke='grey',...){
    isochron_PD(x,nuclide='Rb87',oerr=oerr,sigdig=sigdig,
                show.numbers=show.numbers,levels=levels,
                clabel=clabel,ellipse.fill=ellipse.fill,
                ellipse.stroke=ellipse.stroke,inverse=inverse,
                ci.col=ci.col,line.col=line.col,lwd=lwd,plot=plot,
                exterr=exterr,model=model,wtype=wtype,
                show.ellipses=show.ellipses,hide=hide,omit=omit,
                omit.fill=omit.fill,omit.stroke=omit.stroke,...)
}
#' @rdname isochron
#' @export
isochron.ReOs <- function(x,oerr=3,sigdig=2,show.numbers=FALSE,levels=NA,
                          clabel="",ellipse.fill=c("#00FF0080","#FF000080"),
                          ellipse.stroke='black',inverse=FALSE,
                          ci.col='gray80',line.col='black',lwd=1,
                          plot=TRUE,exterr=TRUE,model=1,wtype='intercept',
                          show.ellipses=1*(model!=2),hide=NULL,
                          omit=NULL,omit.fill=NA,omit.stroke='grey',...){
    isochron_PD(x,nuclide='Re187',oerr=oerr,sigdig=sigdig,
                show.numbers=show.numbers,levels=levels,
                clabel=clabel,ellipse.fill=ellipse.fill,
                ellipse.stroke=ellipse.stroke,inverse=inverse,
                ci.col=ci.col,line.col=line.col,lwd=lwd,plot=plot,
                exterr=exterr,model=model,wtype=wtype,
                show.ellipses=show.ellipses,hide=hide,omit=omit,
                omit.fill=omit.fill,omit.stroke=omit.stroke,...)
}
#' @rdname isochron
#' @export
isochron.SmNd <- function(x,oerr=3,sigdig=2,show.numbers=FALSE,levels=NA,
                          clabel="",ellipse.fill=c("#00FF0080","#FF000080"),
                          ellipse.stroke='black',inverse=FALSE,
                          ci.col='gray80',line.col='black',lwd=1,
                          plot=TRUE,exterr=TRUE,model=1,wtype='intercept',
                          show.ellipses=1*(model!=2),hide=NULL,
                          omit=NULL,omit.fill=NA,omit.stroke='grey',...){
    isochron_PD(x,nuclide='Sm147',oerr=oerr,sigdig=sigdig,
                show.numbers=show.numbers, levels=levels,
                clabel=clabel,ellipse.fill=ellipse.fill,
                ellipse.stroke=ellipse.stroke,inverse=inverse,
                ci.col=ci.col,line.col=line.col,lwd=lwd,plot=plot,
                exterr=exterr,model=model,wtype=wtype,
                show.ellipses=show.ellipses,hide=hide,omit=omit,
                omit.fill=omit.fill,omit.stroke=omit.stroke,...)
}
#' @rdname isochron
#' @export
isochron.LuHf <- function(x,oerr=3,sigdig=2,show.numbers=FALSE,levels=NA,
                          clabel="",ellipse.fill=c("#00FF0080","#FF000080"),
                          ellipse.stroke='black',inverse=FALSE,
                          ci.col='gray80',line.col='black',lwd=1,
                          plot=TRUE,exterr=TRUE,model=1,wtype='intercept',
                          show.ellipses=1*(model!=2),hide=NULL,
                          omit=NULL,omit.fill=NA,omit.stroke='grey',...){
    isochron_PD(x,nuclide='Lu176',oerr=oerr,sigdig=sigdig,
                show.numbers=show.numbers,levels=levels,
                clabel=clabel,ellipse.fill=ellipse.fill,
                ellipse.stroke=ellipse.stroke,inverse=inverse,
                ci.col=ci.col,line.col=line.col,lwd=lwd,plot=plot,
                exterr=exterr,model=model,wtype=wtype,
                show.ellipses=show.ellipses,hide=hide,omit=omit,
                omit.fill=omit.fill,omit.stroke=omit.stroke,...)
}
#' @rdname isochron
#' @export
isochron.UThHe <- function(x,sigdig=2,oerr=3,show.numbers=FALSE,levels=NA,
                           clabel="",ellipse.fill=c("#00FF0080","#FF000080"),
                           ellipse.stroke='black',ci.col='gray80',
                           line.col='black',lwd=1,plot=TRUE,model=1,
                           wtype='intercept',show.ellipses=2*(model!=2),
                           hide=NULL,omit=NULL,omit.fill=NA,
                           omit.stroke='grey',...){
    y <- data2york(x)
    d2calc <- clear(y,hide,omit)
    out <- regression(d2calc,model=model,wtype=wtype)
    out$y0[c('y','s[y]')] <- out$a
    out$age[c('t','s[t]')] <- out$b
    if (inflate(out)){
        out$age['disp[t]'] <- sqrt(out$mswd)*out$age['s[t]']
        out$y0['disp[y]'] <- sqrt(out$mswd)*out$y0['s[y]']
    }
    out$y0label <- quote('He'[0]*'=')
    if (plot){
        scatterplot(y,oerr=oerr,show.numbers=show.numbers,levels=levels,
                    clabel=clabel,ellipse.fill=ellipse.fill,
                    ellipse.stroke=ellipse.stroke,fit=out,
                    show.ellipses=show.ellipses,ci.col=ci.col,
                    line.col=line.col,lwd=lwd,hide=hide,omit=omit,
                    omit.fill=omit.fill,omit.stroke=omit.stroke,...)
        graphics::title(isochrontitle(out,sigdig=sigdig,oerr=oerr,type='U-Th-He'),
                        xlab="P",ylab="He")
    }
    invisible(out)
}
#' @rdname isochron
#' @export
isochron.ThU <- function (x,type=2,oerr=3,sigdig=2,
                          show.numbers=FALSE,levels=NA,clabel="",
                          ellipse.fill=c("#00FF0080","#FF000080"),
                          ellipse.stroke='black',ci.col='gray80',
                          line.col='black',lwd=1,plot=TRUE,
                          exterr=TRUE,model=1,wtype='a',
                          show.ellipses=1*(model!=2),
                          hide=NULL,omit=NULL,omit.fill=NA,
                          omit.stroke='grey',y0option=4,...){
    displabel <- 'dispersion = '
    dispunits <- ''
    if (x$format %in% c(1,2)){
        out <- isochron_ThU_3D(x,type=type,model=model,wtype=wtype,exterr=exterr,
                               hide=hide,omit=omit,y0option=y0option)
        if (model==3){
            if (wtype=='a'){
                displabel <- quote('('^234*'U/'^238*'U)'[a]*'-dispersion = ')
            } else if (wtype=='b'){
                displabel <- quote('('^234*'U/'^232*'Th)-dispersion = ')
            } else if (wtype=='A'){
                displabel <- quote('('^230*'Th/'^238*'U)'[a]*'-dispersion = ')
            } else if (wtype=='B'){
                displabel <- quote('('^230*'Th/'^232*'Th)-dispersion = ')
            }
        }
        dispunits <- ''
    } else if (x$format %in% c(3,4)){
        out <- isochron_ThU_2D(x,type=type,model=model,wtype=wtype,
                               exterr=exterr,hide=hide,omit=omit)
        if (model==3){
            if (type==1 && wtype%in%c('slope',1,'b')){
                age2disp <- TRUE
                Th230U238 <- out$b[1]
            } else if (type==2 && wtype%in%c('intercept',0,'a')){
                age2disp <- TRUE
                Th230U238 <- out$a[1]
            } else {
                age2disp <- FALSE
            }
            if (age2disp){
                l0 <- lambda('Th230')[1]
                dtd08 <- 1/(l0*(1-Th230U238))
                out$disp <- out$disp*dtd08
                dispunits <- ' ka'
            } else {
                displabel <- quote('('^230*'Th/'^232*'Th)'[0]*'-dispersion = ')
            }
        }
    } else {
        stop('Illegal Th-U data format.')
    }
    if (plot){
        scatterplot(out$xyz,oerr=oerr,show.numbers=show.numbers,
                    levels=levels,clabel=clabel,ellipse.fill=ellipse.fill,
                    ellipse.stroke=ellipse.stroke,fit=out,
                    show.ellipses=show.ellipses,ci.col=ci.col,
                    line.col=line.col,lwd=lwd,hide=hide,omit=omit,
                    omit.fill=omit.fill,omit.stroke=omit.stroke,...)
        graphics::title(isochrontitle(out,oerr=oerr,sigdig=sigdig,
                                      type='Th-U',units=' ka',
                                      displabel=displabel,dispunits=dispunits),
                        xlab=out$xlab,ylab=out$ylab)
    }
    invisible(out)
}

isochron_ThU_2D <- function(x,type=2,model=1,wtype='a',
                            exterr=TRUE,hide=NULL,omit=NULL){
    y <- data2york(x,type=type)
    d2calc <- clear(y,hide,omit)
    out <- regression(d2calc,model=model,type="york",wtype=wtype)
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
        get.Th230Th232_0(out$age['t'],Th230Th232[1],Th230Th232[2])
    if (inflate(out)){
        out$age['disp[t]'] <- get.ThU.age(Th230U238[1],
                                          sqrt(out$mswd)*Th230U238[2],
                                          exterr=exterr)['s[t]']
        out$y0['disp[y]'] <- get.Th230Th232_0(out$age['t'],
                                              Th230Th232[1],
                                              sqrt(out$mswd)*Th230Th232[2])[2]
    }
    out$xlab <- x.lab
    out$ylab <- y.lab
    out$y0label <- quote('('^230*'Th/'^232*'Th)'[0]*'=')
    out
}
isochron_ThU_3D <- function(x,type=2,model=1,wtype='a',exterr=TRUE,
                            hide=NULL,omit=NULL,y0option=4){
    if (type == 1){ # 0/2 vs 8/2
        osmond <- FALSE
        ia <- 'A'
        ib <- 'B'
        i48 <- 'b'
        i08 <- 'B'
        i02 <- 'A'
        id <- c('X','sX','Z','sZ','rXZ')
        x.lab <- quote(''^238*'U/'^232*'Th')
        y.lab <- quote(''^230*'Th/'^232*'Th')
    } else if (type == 2){ # 0/8 vs 2/8
        osmond <- TRUE
        ia <- 'A'
        ib <- 'B'
        i48 <- 'a'
        i08 <- 'A'
        i02 <- 'B'
        id <- c('X','sX','Z','sZ','rXZ')
        x.lab <- quote(''^232*'Th/'^238*'U')
        y.lab <- quote(''^230*'Th/'^238*'U')
    } else if (type == 3){ # 4/2 vs 8/2
        osmond <- FALSE
        ia <- 'a'
        ib <- 'b'
        i48 <- 'b'
        i08 <- 'B'
        i02 <- 'A'
        id <- c('X','sX','Y','sY','rXY')
        x.lab <- quote(''^238*'U/'^232*'Th')
        y.lab <- quote(''^234*'U/'^232*'Th')
    } else if (type == 4){ # 4/8 vs 2/8
        osmond <- TRUE
        ia <- 'a'
        ib <- 'b'
        i48 <- 'a'
        i08 <- 'A'
        i02 <- 'B'
        id <- c('X','sX','Y','sY','rXY')
        x.lab <- quote(''^232*'Th/'^238*'U')
        y.lab <- quote(''^234*'U/'^238*'U')
    }
    tit <- data2tit(x,osmond=osmond)
    d2calc <- clear(tit,hide,omit)
    out <- regression(d2calc,model=model,type="titterington",wtype=wtype)
    out$xyz <- tit
    out$a <- c(out$par[ia],sqrt(out$cov[ia,ia]))
    out$b <- c(out$par[ib],sqrt(out$cov[ib,ib]))
    out$cov.ab <- out$cov[ia,ib]
    out$PAR <- out$par[c(i48,i08,i02)]
    out$COV <- out$cov[c(i48,i08,i02),c(i48,i08,i02)]
    names(out$PAR) <- rownames(out$COV) <- colnames(out$COV) <- c('i48','i08','i02')
    tst <- get.ThU.age(out$par[i08],sqrt(out$cov[i08,i08]),
                       out$par[i48],sqrt(out$cov[i48,i48]),
                       out$cov[i48,i08],exterr=exterr,jacobian=TRUE)
    out$age['t'] <- tst['t']
    out$age['s[t]'] <- tst['s[t]']
    if (inflate(out)){
        tst <- get.ThU.age(out$par[i08],sqrt(out$mswd*out$cov[i08,i08]),
                           out$par[i48],sqrt(out$mswd*out$cov[i48,i48]),
                           out$mswd*out$cov[i48,i08],exterr=exterr,jacobian=TRUE)
        out$age['disp[t]'] <- tst['s[t]']
    }
    out <- getThUy0(out=out,tst=tst,option=y0option,exterr=exterr)
    out$xlab <- x.lab
    out$ylab <- y.lab
    out$xyz <- subset(out$xyz,select=id)
    out
}
getThUy0 <- function(out,tst,option=1,exterr=FALSE){
    if (option==1){
        out$y0['y'] <- out$PAR['i48']
        out$y0['s[y]'] <- sqrt(out$COV['i48','i48'])
        out$y0label <- quote('('^234*'U/'^238*'U)'[a]*'=')
    } else if (option==2){
        out$y0['y'] <- out$PAR['i02']
        out$y0['s[y]'] <- sqrt(out$COV['i02','i02'])
        out$y0label <- quote('('^230*'Th/'^232*'Th)'[d]*'=')
    } else if (option==3){
        out$y0['y'] <- out$PAR['i08']
        out$y0['s[y]'] <- sqrt(out$COV['i08','i08'])
        out$y0label <- quote('('^230*'Th/'^238*'U)'[a]*'=')
    } else {
        l4 <- lambda('U234')
        out$y0['y'] <- 1 + (out$PAR['i48']-1)*exp(l4[1]*tst['t'])
        E <- matrix(0,2,2)
        E[1,1] <- out$COV['i48','i48']
        E[2,2] <- l4[2]^2
        J <- matrix(0,3,2)
        J[1,1] <- tst['dt.d48']
        J[2,1] <- 1
        J[3,2] <- ifelse(exterr,1,0)
        E2 <- J %*% E %*% t(J)
        J2 <- rep(0,3)
        J2[1] <- l4[1]*(out$PAR['i48']-1)*exp(l4[1]*tst['t'])
        J2[2] <- exp(l4[1]*tst['t'])
        J2[3] <- ifelse(exterr,(out$PAR['i48']-1)*exp(l4[1]*tst['t'])*tst['t'],0)
        out$y0['s[y]'] <- sqrt(J2 %*% E2 %*% J2)
        out$y0label <- quote('('^234*'U/'^238*'U)'[i]*'=')
    }
    if (inflate(out)){ # overwrite dispersion to remove decay constant errors
        if (option<4){
            out$age['disp[t]'] <- sqrt(out$mswd)*out$age['s[t]']
            out$age['disp[y]'] <- sqrt(out$mswd)*out$age['s[y]']
        } else {
            E[1,1] <- out$mswd*out$COV['i48','i48']
            E3 <- J2 %*% (J %*% E %*% t(J)) %*% J2
            out$age['disp[t]'] <- sqrt(E3[1,1])
            out$y0['disp[y]'] <- sqrt(E3[2,2])
        }
    }
    out
}

isochron_PD <- function(x,nuclide,oerr=3,sigdig=2,
                        show.numbers=FALSE,levels=NA,clabel="",
                        ellipse.fill=c("#00FF0080","#FF000080"),
                        ellipse.stroke='black',inverse=FALSE,
                        ci.col='gray80',line.col='black',lwd=1,
                        plot=TRUE,exterr=TRUE,model=1,wtype=NA,
                        show.ellipses=1*(model!=2),bratio=1,
                        hide=NULL,omit=NULL,...){
    y <- data2york(x,inverse=inverse)
    d2calc <- clear(y,hide,omit)
    out <- regression(d2calc,model=model,wtype=wtype)
    out$y0[c('y','s[y]')] <- out$a
    if (inverse){
        DP <- -out$b[1]/out$a[1]
        sDP <- DP*sqrt((out$a[2]/out$a[1])^2 + (out$b[2]/out$b[1])^2 -
                       2*out$cov.ab/(out$a[1]*out$b[1]))
    } else {
        DP <- out$b[1]
        sDP <- out$b[2]
    }
    out$age[c('t','s[t]')] <- get.PD.age(DP,sDP,nuclide,
                                         exterr=exterr,bratio=bratio)
    if (inflate(out)){
        out$age['disp[t]'] <- get.PD.age(DP,sqrt(out$mswd)*sDP,nuclide,
                                         exterr=exterr,bratio=bratio)[2]
        out$y0['disp[y]'] <- sqrt(out$mswd)*out$y0['s[y]']
    }
    if (wtype%in%c('slope',1,'b')){
        dDPdt <- lambda(nuclide)[1]*(1+DP)
        dDPdb <- ifelse(inverse,1/out$a[1],1)
        out$disp <- out$disp*dDPdb/dDPdt
        dispunits <- ' Ma'
    } else {
        dispunits <- ''
    }
    lab <- get.isochron.labels(nuclide=nuclide,inverse=inverse)
    out$y0label <- substitute(a*b*c,list(a='(',b=lab$y,c=quote(')'[0]*'=')))
    if (plot){
        scatterplot(y,oerr=oerr,show.ellipses=show.ellipses,
                    show.numbers=show.numbers,levels=levels,
                    clabel=clabel,ellipse.fill=ellipse.fill,
                    ellipse.stroke=ellipse.stroke,fit=out,
                    ci.col=ci.col,line.col=line.col,lwd=lwd,
                    hide=hide,omit=omit,...)
        graphics::title(isochrontitle(out,oerr=oerr,sigdig=sigdig,
                                      type='PD',dispunits=dispunits),
                        xlab=lab$x,ylab=lab$y)
    }
    invisible(out)
}

get.isochron.labels <- function(nuclide,inverse=FALSE){
    out <- list()
    if (identical(nuclide,'Th232')){
        if (inverse){
            out$x <- quote(''^232*'Th/'^208*'Pb')
            out$y <- quote(''^204*'Pb/'^208*'Pb')
        } else {
            out$x <- quote(''^232*'Th/'^204*'Pb')
            out$y <- quote(''^208*'Pb/'^204*'Pb')
        }
    } else if (identical(nuclide,'Sm147')){
        if (inverse){
            out$x <- quote(''^147*'Sm/'^143*'Nd')
            out$y <- quote(''^144*'Nd/'^143*'Nd')
        } else {
            out$x <- quote(''^147*'Sm/'^144*'Nd')
            out$y <- quote(''^143*'Nd/'^144*'Nd')
        }
    } else if (identical(nuclide,'Re187')){
        if (inverse){
            out$x <- quote(''^187*'Re/'^187*'Os')
            out$y <- quote(''^188*'Os/'^187*'Os')
        } else {
            out$x <- quote(''^187*'Re/'^188*'Os')
            out$y <- quote(''^187*'Os/'^188*'Os')
        }
    } else if (identical(nuclide,'Rb87')){
        if (inverse){
            out$x <- quote(''^87*'Rb/'^87*'Sr')
            out$y <- quote(''^86*'Sr/'^87*'Sr')
        } else {
            out$x <- quote(''^87*'Rb/'^86*'Sr')
            out$y <- quote(''^87*'Sr/'^86*'Sr')
        }
    } else if (identical(nuclide,'Lu176')){
        if (inverse){
            out$x <- quote(''^176*'Lu/'^176*'Hf')
            out$y <- quote(''^177*'Hf/'^176*'Hf')
        } else {
            out$x <- quote(''^176*'Lu/'^177*'Hf')
            out$y <- quote(''^176*'Hf/'^177*'Hf')
        }
    } else if (identical(nuclide,'K40')){
        if (inverse){
            out$x <- quote(''^40*'K/'^40*'Ca')
            out$y <- quote(''^44*'Ca/'^40*'Ca')
        } else {
            out$x <- quote(''^40*'K/'^44*'Ca')
            out$y <- quote(''^40*'Ca/'^44*'Ca')
        }
    }
    out
}

get.limits <- function(x,sx){
    minx <- min(x-3*sx,na.rm=TRUE)
    maxx <- max(x+3*sx,na.rm=TRUE)
    c(minx,maxx)
}

plot_PbPb_evolution <- function(from=0,to=4570,inverse=TRUE){
    nn <- 50
    tijd <- seq(from=from,to=to,length.out=nn)
    ticks <- pretty(tijd)
    tijd[nn] <- max(ticks)
    xy <- stacey.kramers(tijd,inverse=inverse)
    graphics::lines(xy[,1],xy[,2])
    xy <- stacey.kramers(ticks,inverse=inverse)
    graphics::points(xy[,1],xy[,2],pch=20)
    graphics::text(xy[,1],xy[,2],labels=ticks,pos=3)
}

isochrontitle <- function(fit,oerr=3,sigdig=2,type=NULL,
                          units=' Ma',displabel='dispersion =',
                          dispunits='',ski=NULL,y0option=1,...){
    content <- list()
    if (is.null(type)){
        content[[1]] <- maintit(x=fit$a[1],sx=fit$a[-1],n=fit$n,
                                units=units,prefix='intercept =',
                                sigdig=sigdig,oerr=oerr,df=fit$df)
        content[[2]] <- maintit(x=fit$b[1],sx=fit$b[-1],ntit='',
                                units=units,prefix='slope =',
                                sigdig=sigdig,oerr=oerr,df=fit$df)
    } else if (type=='U-Pb'){
        if (is.null(fit$posterior) || 't'%ni%names(fit$posterior)){
            content[[1]] <- maintit(x=fit$age[1],sx=fit$age[-1],n=fit$n,
                                    units=units,sigdig=sigdig,
                                    oerr=oerr,df=fit$df)
        } else {
            content[[1]] <- bayestit(x=fit$par['t'],XL=fit$posterior$t,
                                     n=fit$n,sigdig=sigdig,oerr=oerr)
        }
        if(is.null(fit$posterior)) pnames <- NULL
        else pnames <- names(fit$posterior)
        if (is.null(pnames)) ipar <- NULL
        else if (y0option==2 && 'U48i'%in%pnames) ipar <- 'U48i'
        else if (y0option==3 && 'ThUi'%in%pnames) ipar <-'ThUi'
        else ipar <- NULL
        if (is.null(ipar)){
            content[[2]] <- maintit(x=fit$y0[1],sx=fit$y0[-1],ntit='',
                                    units='',prefix=fit$y0label,
                                    sigdig=sigdig,oerr=oerr,df=fit$df)
        } else {
            content[[2]] <- bayestit(x=fit$par[ipar],XL=fit$posterior[[ipar]],
                                     ntit='',sigdig=sigdig,oerr=oerr,units='',
                                     prefix=fit$y0label)
        }
    } else {
        content[[1]] <- maintit(x=fit$age[1],sx=fit$age[-1],n=fit$n,
                                units=units,sigdig=sigdig,
                                oerr=oerr,df=fit$df)
        content[[2]] <- maintit(x=fit$y0[1],sx=fit$y0[-1],ntit='',
                                units='',prefix=fit$y0label,
                                sigdig=sigdig,oerr=oerr,df=fit$df)
    }
    if (fit$model==1){
        content[[3]] <- mswdtit(mswd=fit$mswd,p=fit$p.value,sigdig=sigdig)
    } else if (fit$model%in%c(3,5)){
        content[[3]] <- disptit(w=fit$disp[1],sw=fit$disp[2],sigdig=sigdig,
                                oerr=oerr,prefix=displabel,units=dispunits)
    }
    nl <- length(content)
    if (!is.null(ski)){
        growthline <- paste0('intercepts growth curve at ',
                             roundit(ski[1],sigdig=sigdig,text=TRUE))
        if (length(ski)>1){
            growthline <- paste0(growthline,' and ',
                                 roundit(ski[2],sigdig=sigdig,text=TRUE))
        }
        growthline <- paste0(growthline,' Ma')
        nl <- nl + 1
        content[[nl]] <- growthline
    }
    for (i in nl:1){
        mymtext(content[[i]],line=nl-i,...)
    }
}
