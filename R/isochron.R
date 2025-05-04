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
#' If \code{x} has class \code{PbPb}, \code{ArAr} or \code{PD},
#' \code{wtype} can have one of two values:
#'
#' \itemize{
#' \item \code{1}: attribute the overdispersion to variability in the
#' non-radiogenic component, as controlled by \code{settings('iratio',...)}
#'
#' \item \code{2}: attribute the overdispersion to variability in the
#' age, i.e. to diachronous closure of the isotope system.
#' }
#'
#' otherwise, \code{wtype} can have one of four values:
#'
#' \itemize{
#'
#' \item \code{1}: attributes the overdispersion to the y-intercept of
#' the equivalent conventional isochron.
#'
#' \item \code{2}: attributes the overdispersion to the slope of the
#' equivalent conventional isochron.
#'
#' \item \code{'A'}: only available if \code{x} has class \code{ThU} and
#' \code{x$format} is 1 or 2. Attributes the overdispersion to the
#' authigenic \eqn{^{230}}Th/\eqn{^{238}}U-intercept of the isochron.
#'
#' \item \code{'B'}: only available if \code{x} has class \code{ThU} and
#' \code{x$format} is 1 or 2. Attributes the overdispersion to the
#' \eqn{^{230}}Th/\eqn{^{232}}Th-slope of the isochron.
#' }
#'
#' @param anchor control parameters to fix the intercept age or
#'     non-radiogenic composition of the isochron fit. This can be a
#'     scalar or a vector.
#'
#' If \code{anchor[1]=0}: do not anchor the isochron.
#'
#' If \code{anchor[1]=1}: fix the non-radiogenic composition at the
#' values stored in \code{settings('iratio',...)}, OR, if \code{x} has
#' class \code{other}, fix the intercept at the value stored in
#' \code{anchor[2]}.
#'
#' If \code{anchor[1]=2}: fix the age at the value stored in \code{anchor[2]}.
#'
#' If \code{x} has class \code{UPb} and \code{anchor[1]=3}: anchor the
#' non-radiogenic component to the Stacey-Kramers mantle evolution
#' model.
#'
#' @param flippable controls if generic data (where \code{x} has class
#'     \code{other} and \code{x$format} is either \code{4} or
#'     \code{5}) should be treated as inverse isochrons
#'     (\code{flippable=1}) or as conventional isochrons
#'     (\code{flippable=2}). If \code{flippable=0} (which is the
#'     default value), then the data are passed on to
#'     \code{isochron.default}.
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
#' @param ... optional arguments to be passed on to \code{\link{scatterplot}}
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
#' \item{df}{the degrees of freedom of the linear fit (\eqn{df=n-2}
#' for non-anchored fits)}
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
#' the standard deviation of the (assumed) normally distributed
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
#' OR, if \code{x} has class \code{other} and \code{x$format} is
#' either \code{4} or \code{5} and \code{flippable} is not \code{0},
#' returns
#'
#' \code{Dd}: the ratio of the inherited radiogenic daughter to its
#' nonradiogenic sister isotope
#'
#' \code{DP}: the ratio fo the radiogic daughter to its radioactive parent
#'
#' \code{cov.DdDP}: the covariance between \code{Dd} and \code{DP}.
#'
#' In the remaining types of \code{other} data, the intercept \code{a}
#' and \code{b} are returned along with their covariance.
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
                             levels=NULL,clabel="",xlab='x',ylab='y',
                             ellipse.fill=c("#00FF0080","#FF000080"),
                             ellipse.stroke='black',ci.col='gray80',
                             line.col='black',lwd=1,plot=TRUE,title=TRUE,
                             model=1,wtype=1,anchor=0,show.ellipses=1*(model!=2),
                             hide=NULL,omit=NULL,omit.fill=NA,
                             omit.stroke='grey',...){
    d2calc <- clear(x,hide,omit)
    if (model>1 | anchor[1]==2){
        fit <- MLyork(d2calc,anchor=anchor,model=model,wtype=wtype)
        if (model==3) fit$disp <- fit$w
    } else {
        if (anchor[1]<1){
            fit <- regression(d2calc)
        } else {
            if (length(anchor>1)) y0 <- anchor[2]
            else stop("anchor must be a vector of at least two numbers.")
            sy0 <- ifelse(length(anchor)>2,anchor[3],0)
            fit <- anchoredYork(d2calc,y0=y0,sy0=sy0)
        }
    }
    out <- ab2y0t(x=x,fit=fit)
    genericisochronplot(x=x,fit=out,oerr=oerr,sigdig=sigdig,
                        show.numbers=show.numbers,levels=levels,clabel=clabel,
                        xlab=xlab,ylab=ylab,ellipse.fill=ellipse.fill,
                        ellipse.stroke=ellipse.stroke,ci.col=ci.col,
                        line.col=line.col,lwd=lwd,plot=plot,title=title,
                        show.ellipses=show.ellipses,hide=hide,omit=omit,
                        omit.fill=omit.fill,omit.stroke=omit.stroke,...)
}
#' @rdname isochron
#' @export
isochron.other <- function(x,oerr=3,sigdig=2,show.numbers=FALSE,
                           levels=NULL,clabel="",xlab='x',ylab='y',
                           ellipse.fill=c("#00FF0080","#FF000080"),
                           ellipse.stroke='black',ci.col='gray80',
                           line.col='black',lwd=1,plot=TRUE,
                           title=TRUE, model=1,wtype=1,anchor=0,
                           flippable=0,show.ellipses=1*(model!=2),
                           hide=NULL,omit=NULL,omit.fill=NA,
                           omit.stroke='grey',...){
    if (x$format<4){
        stop('Invalid format for isochron regression.')
    } else if (x$format==6){
        d2calc <- clear(x,hide,omit)
        fit <- regression(d2calc$x,model=model,type='ogls')
        out <- ab2y0t(x=d2calc$x,fit=fit)
        genericisochronplot(x=x,fit=out,oerr=oerr,sigdig=sigdig,
                            show.numbers=show.numbers,levels=levels,clabel=clabel,
                            xlab=xlab,ylab=ylab,ellipse.fill=ellipse.fill,
                            ellipse.stroke=ellipse.stroke,ci.col=ci.col,
                            line.col=line.col,lwd=lwd,plot=plot,title=title,
                            show.ellipses=show.ellipses,hide=hide,omit=omit,
                            omit.fill=omit.fill,omit.stroke=omit.stroke,...)
    } else if (flippable%in%c(1,2)){ # treat as geochronology data
        wtype <- checkWtype(wtype=wtype,anchor=anchor,model=model)
        inverse <- (flippable==1)
        fit <- flipper(x,inverse=inverse,model=model,wtype=wtype,
                       anchor=anchor,hide=hide,omit=omit)
        out <- ab2y0t(x=x,fit=fit,inverse=inverse,wtype=wtype)
        scatterplot(out$xyz,oerr=oerr,show.ellipses=show.ellipses,
                    show.numbers=show.numbers,levels=levels,
                    clabel=clabel,ellipse.fill=ellipse.fill,
                    ellipse.stroke=ellipse.stroke,fit=out,
                    ci.col=ci.col,line.col=line.col,lwd=lwd,
                    hide=hide,omit=omit,omit.fill=omit.fill,
                    omit.stroke=omit.stroke,...)
        showDispersion(out,inverse=(flippable==1),wtype=wtype)
        if (title){
            main <- isochrontitle(out,oerr=oerr,sigdig=sigdig,
                                  units='',type='generic')
        } else {
            main <- NULL
        }
        graphics::title(main=main,xlab=xlab,ylab=ylab)
    } else { # general purpose regression
        yd <- data2york(x)
        out <- isochron(yd,oerr=oerr,sigdig=sigdig,
                        show.numbers=show.numbers,levels=levels,
                        clabel=clabel,xlab=xlab,ylab=ylab,
                        ellipse.fill=ellipse.fill,
                        ellipse.stroke=ellipse.stroke,ci.col=ci.col,
                        line.col=line.col,lwd=lwd,plot=plot,
                        title=title,model=model,wtype=wtype,anchor=anchor,
                        show.ellipses=show.ellipses,hide=hide,omit=omit,
                        omit.fill=omit.fill,omit.stroke=omit.stroke,...)
    }
    invisible(out)
}
genericisochronplot <- function(x,fit,oerr=3,sigdig=2,show.numbers=FALSE,
                                levels=NULL,clabel="",xlab='x',ylab='y',
                                ellipse.fill=c("#00FF0080","#FF000080"),
                                ellipse.stroke='black',ci.col='gray80',
                                line.col='black',lwd=1,plot=TRUE,
                                title=TRUE,show.ellipses=TRUE,
                                hide=NULL,omit=NULL,omit.fill=NA,
                                omit.stroke='grey',...){
    if (plot){
        y <- data2york(x)
        scatterplot(y,oerr=oerr,show.ellipses=show.ellipses,
                    show.numbers=show.numbers,levels=levels,
                    clabel=clabel,ellipse.fill=ellipse.fill,
                    ellipse.stroke=ellipse.stroke,fit=fit,
                    ci.col=ci.col,line.col=line.col,lwd=lwd,
                    hide=hide,omit=omit,omit.fill=omit.fill,
                    omit.stroke=omit.stroke,...)
        if (title){
            main <- isochrontitle(fit,oerr=oerr,sigdig=sigdig,units='')
        } else {
            main <- NULL
        }
        graphics::title(main=main,xlab=xlab,ylab=ylab)
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
#' @param taxis logical. If \code{TRUE}, replaces the x-axis of the
#'     inverse isochron with a time scale. Only applies if
#'     \code{inverse} is \code{TRUE} or, when \code{x} has class
#'     \code{U-Pb}, if \code{x$format} is 4 or higher.
#'
#' @rdname isochron
#' @export
isochron.UPb <- function(x,oerr=3,sigdig=2,show.numbers=FALSE,
                         levels=NULL,clabel="",joint=TRUE,
                         ellipse.fill=c("#00FF0080","#FF000080"),
                         ellipse.stroke='black',type=1,
                         ci.col='gray80',line.col='black',lwd=1,
                         plot=TRUE,title=TRUE,exterr=FALSE,model=1,
                         show.ellipses=1*(model!=2),anchor=0,
                         hide=NULL,omit=NULL,omit.fill=NA,
                         omit.stroke='grey',y0option=1,taxis=FALSE,...){
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
        type <- checkIsochronType(x,type=type)
        ns <- length(x)
        calcit <- (1:ns)%ni%c(hide,omit)
        x2calc <- subset(x,subset=calcit)
        fit <- ludwig(x2calc,model=model,anchor=anchor,
                      exterr=exterr,type=ifelse(joint,0,type))
        out <- ab2y0t(x=x,fit=fit,type=type,exter=exterr,y0option=y0option)
        if (plot){
            scatterplot(out$XY,oerr=oerr,show.ellipses=show.ellipses,
                        show.numbers=show.numbers,levels=levels,
                        clabel=clabel,ellipse.fill=ellipse.fill,
                        ellipse.stroke=ellipse.stroke,fit=out,
                        ci.col=ci.col,line.col=line.col,lwd=lwd,
                        hide=hide,omit=omit,omit.fill=omit.fill,
                        omit.stroke=omit.stroke,taxis=taxis,...)
            if (taxis) add_taxis(x=x,fit=out,type=type)
            if (!joint | x$format>8){
                showDispersion(out,inverse=TRUE,wtype=anchor[1],type=type)
            }
            dispunits <- getDispUnits_UPb(x=x,joint=joint,anchor=anchor)
            if (title){
                main <- isochrontitle(out,oerr=oerr,sigdig=sigdig,type='U-Pb',
                                      y0option=y0option,dispunits=dispunits)
            } else {
                main <- NULL
            }
            lab <- getIsochronLabels(x=x,type=type,taxis=taxis)
            graphics::title(main=main,xlab=lab$x,ylab=lab$y)
        }
    }
    invisible(out)
}

checkIsochronType <- function(x,type=1){
    if (x$format%in%c(9,11,119) & type%ni%c(1,3)) return(1)
    else if (x$format%in%c(10,12,1210) & type%ni%c(2,4)) return(2)
    else return(type)
}
checkWtype <- function(wtype=1,anchor=0,model=1){
    if (anchor[1]==1 & model==3) return(1)
    else if (anchor[1]==2 & model==3) return(2)
    else return(wtype)
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
isochron.PbPb <- function(x,oerr=3,sigdig=2,show.numbers=FALSE,levels=NULL,
                          clabel="",ellipse.fill=c("#00FF0080","#FF000080"),
                          ellipse.stroke='black',inverse=TRUE,
                          ci.col='gray80',line.col='black',lwd=1,
                          plot=TRUE,title=TRUE,exterr=FALSE,model=1,
                          wtype=1,anchor=0,growth=FALSE,show.ellipses=1*(model!=2),
                          hide=NULL,omit=NULL,omit.fill=NA,omit.stroke='grey',...){
    wtype <- checkWtype(wtype=wtype,anchor=anchor,model=model)
    fit <- flipper(x,inverse=inverse,model=model,wtype=wtype,
                   anchor=anchor,hide=hide,omit=omit,type="d")
    out <- ab2y0t(x=x,fit=fit,inverse=inverse,exterr=exterr,wtype=wtype)
    dispunits <- getDispUnits(model=model,wtype=wtype,anchor=anchor)
    if (plot) {
        scatterplot(out$xyz,oerr=oerr,show.ellipses=show.ellipses,
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
            out$ski <- SK_intersection(out,inverse=inverse)
        } else {
            out$ski <- NULL
        }
        showDispersion(out,inverse=inverse,wtype=wtype,type='d')
        if (title){
            main <- isochrontitle(out,oerr=oerr,sigdig=sigdig,type='Pb-Pb',
                                  dispunits=dispunits,ski=out$ski)
        } else {
            main <- NULL
        }
        lab <- getIsochronLabels(x=x,inverse=inverse)
        graphics::title(main=main,xlab=lab$x,ylab=lab$y)
    }
    invisible(out)
}
SK_intersection <- function(fit,inverse,m=0,M=5000){
    SKi_misfit <- function(tt,a,b){
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
    } else if (sign(SKi_misfit(m,a,b)) == sign(SKi_misfit(M,a,b))){
        ski1 <- SK_intersection(fit,inverse,m=m,M=m+(M-m)/2)
        ski2 <- SK_intersection(fit,inverse,m=m+(M-m)/2,M=M)
        out <- unique(c(ski1,ski2))
    } else {
        out <- stats::uniroot(SKi_misfit,interval=c(m,M),a=a,b=b)$root
    }
    out
}
#' @rdname isochron
#' @export
isochron.ArAr <- function(x,oerr=3,sigdig=2,show.numbers=FALSE,
                          levels=NULL,clabel="",
                          ellipse.fill=c("#00FF0080","#FF000080"),
                          ellipse.stroke='black',inverse=TRUE,
                          ci.col='gray80',line.col='black',lwd=1,
                          plot=TRUE,title=TRUE,exterr=FALSE,model=1,
                          wtype=1,anchor=0,show.ellipses=1*(model!=2),
                          hide=NULL,omit=NULL,omit.fill=NA,
                          omit.stroke='grey',taxis=FALSE,...){
    taxis <- taxis & inverse
    wtype <- checkWtype(wtype=wtype,anchor=anchor,model=model)
    fit <- flipper(x,inverse=inverse,model=model,wtype=wtype,
                   anchor=anchor,hide=hide,omit=omit,type="p",
                   J=x$J[1],sJ=x$J[2])
    out <- ab2y0t(x=x,fit=fit,inverse=inverse,exterr=exterr,wtype=wtype)
    dispunits <- getDispUnits(model=model,wtype=wtype,anchor=anchor)
    if (plot) {
        scatterplot(out$xyz,oerr=oerr,show.ellipses=show.ellipses,
                    show.numbers=show.numbers,levels=levels,
                    clabel=clabel,ellipse.fill=ellipse.fill,
                    ellipse.stroke=ellipse.stroke,fit=out,
                    ci.col=ci.col,line.col=line.col,lwd=lwd,
                    hide=hide,omit=omit,omit.fill=omit.fill,
                    omit.stroke=omit.stroke,taxis=taxis,...)
        if (taxis) add_taxis(x=x,fit=out)
        showDispersion(out,inverse=inverse,wtype=wtype)
        if (title){
            main <- isochrontitle(out,oerr=oerr,sigdig=sigdig,
                                  dispunits=dispunits,type='Ar-Ar')
        } else {
            main <- NULL
        }
        lab <- getIsochronLabels(x=x,inverse=inverse,taxis=taxis)
        graphics::title(main=main,xlab=lab$x,ylab=lab$y)
    }
    invisible(out)
}
#' @rdname isochron
#' @export
isochron.ThPb <- function(x,oerr=3,sigdig=2,show.numbers=FALSE,levels=NULL,
                          clabel="",ellipse.fill=c("#00FF0080","#FF000080"),
                          ellipse.stroke='black',inverse=FALSE,
                          ci.col='gray80',line.col='black',lwd=1,
                          plot=TRUE,title=TRUE,exterr=FALSE,model=1,
                          wtype=1,anchor=0,show.ellipses=1*(model!=2),
                          hide=NULL,omit=NULL,omit.fill=NA,
                          omit.stroke='grey',taxis=FALSE,...){
    isochron_PD(x,oerr=oerr,sigdig=sigdig,show.numbers=show.numbers,
                levels=levels,clabel=clabel,ellipse.fill=ellipse.fill,
                ellipse.stroke=ellipse.stroke,inverse=inverse,
                ci.col=ci.col,line.col=line.col,lwd=lwd,plot=plot,
                title=title,exterr=exterr,model=model,wtype=wtype,
                anchor=anchor,show.ellipses=show.ellipses,hide=hide,
                omit=omit,omit.fill=omit.fill,omit.stroke=omit.stroke,
                taxis=taxis,...)
}
#' @param bratio the \eqn{^{40}}K branching ratio.
#' @rdname isochron
#' @export
isochron.KCa <- function(x,oerr=3,sigdig=2,show.numbers=FALSE,levels=NULL,
                         clabel="",inverse=FALSE,ci.col='gray80',
                         ellipse.fill=c("#00FF0080","#FF000080"),
                         ellipse.stroke='black',line.col='black',
                         lwd=1,plot=TRUE,title=TRUE,exterr=FALSE,model=1,
                         wtype=1,anchor=0,show.ellipses=1*(model!=2),
                         hide=NULL,omit=NULL,omit.fill=NA,omit.stroke='grey',
                         taxis=FALSE,bratio=0.895,...){
    isochron_PD(x,oerr=oerr,sigdig=sigdig,show.numbers=show.numbers,
                levels=levels,clabel=clabel,ellipse.fill=ellipse.fill,
                ellipse.stroke=ellipse.stroke,inverse=inverse,
                ci.col=ci.col,line.col=line.col,lwd=lwd,plot=plot,
                title=title,exterr=exterr,model=model,wtype=wtype,
                anchor=anchor,show.ellipses=show.ellipses,bratio=bratio,
                hide=hide,omit=omit,omit.fill=omit.fill,
                omit.stroke=omit.stroke,taxis=taxis,...)
}
#' @rdname isochron
#' @export
isochron.RbSr <- function(x,oerr=3,sigdig=2,show.numbers=FALSE,levels=NULL,
                          clabel="",ellipse.fill=c("#00FF0080","#FF000080"),
                          ellipse.stroke='black',inverse=FALSE,
                          ci.col='gray80',line.col='black',lwd=1,
                          plot=TRUE,title=TRUE,exterr=FALSE,model=1,wtype=1,
                          anchor=0,show.ellipses=1*(model!=2),hide=NULL,
                          omit=NULL,omit.fill=NA,omit.stroke='grey',
                          taxis=FALSE,...){
    isochron_PD(x,oerr=oerr,sigdig=sigdig,show.numbers=show.numbers,
                levels=levels,clabel=clabel,ellipse.fill=ellipse.fill,
                ellipse.stroke=ellipse.stroke,inverse=inverse,
                ci.col=ci.col,line.col=line.col,lwd=lwd,plot=plot,
                title=title,exterr=exterr,model=model,wtype=wtype,
                anchor=anchor,show.ellipses=show.ellipses,hide=hide,omit=omit,
                omit.fill=omit.fill,omit.stroke=omit.stroke,taxis=taxis,...)
}
#' @rdname isochron
#' @export
isochron.ReOs <- function(x,oerr=3,sigdig=2,show.numbers=FALSE,levels=NULL,
                          clabel="",ellipse.fill=c("#00FF0080","#FF000080"),
                          ellipse.stroke='black',inverse=FALSE,
                          ci.col='gray80',line.col='black',lwd=1,
                          plot=TRUE,title=TRUE,exterr=FALSE,model=1,wtype=1,
                          anchor=0,show.ellipses=1*(model!=2),hide=NULL,omit=NULL,
                          omit.fill=NA,omit.stroke='grey',taxis=FALSE,...){
    isochron_PD(x,oerr=oerr,sigdig=sigdig,show.numbers=show.numbers,
                levels=levels,clabel=clabel,ellipse.fill=ellipse.fill,
                ellipse.stroke=ellipse.stroke,inverse=inverse,
                ci.col=ci.col,line.col=line.col,lwd=lwd,plot=plot,
                title=title,exterr=exterr,model=model,wtype=wtype,
                anchor=anchor,show.ellipses=show.ellipses,hide=hide,omit=omit,
                omit.fill=omit.fill,omit.stroke=omit.stroke,taxis=taxis,...)
}
#' @rdname isochron
#' @export
isochron.SmNd <- function(x,oerr=3,sigdig=2,show.numbers=FALSE,levels=NULL,
                          clabel="",ellipse.fill=c("#00FF0080","#FF000080"),
                          ellipse.stroke='black',inverse=FALSE,
                          ci.col='gray80',line.col='black',lwd=1,
                          plot=TRUE,title=TRUE,exterr=FALSE,model=1,wtype=1,
                          anchor=0,show.ellipses=1*(model!=2),hide=NULL,omit=NULL,
                          omit.fill=NA,omit.stroke='grey',taxis=FALSE,...){
    isochron_PD(x,oerr=oerr,sigdig=sigdig,show.numbers=show.numbers,
                levels=levels,clabel=clabel,ellipse.fill=ellipse.fill,
                ellipse.stroke=ellipse.stroke,inverse=inverse,
                ci.col=ci.col,line.col=line.col,lwd=lwd,plot=plot,
                title=title,exterr=exterr,model=model,wtype=wtype,
                anchor=anchor,show.ellipses=show.ellipses,hide=hide,omit=omit,
                omit.fill=omit.fill,omit.stroke=omit.stroke,taxis=taxis,...)
}
#' @rdname isochron
#' @export
isochron.LuHf <- function(x,oerr=3,sigdig=2,show.numbers=FALSE,levels=NULL,
                          clabel="",ellipse.fill=c("#00FF0080","#FF000080"),
                          ellipse.stroke='black',inverse=FALSE,
                          ci.col='gray80',line.col='black',lwd=1,plot=TRUE,
                          title=TRUE,exterr=FALSE,model=1,wtype=1,anchor=0,
                          show.ellipses=1*(model!=2),hide=NULL,omit=NULL,
                          omit.fill=NA,omit.stroke='grey',taxis=FALSE,...){
    isochron_PD(x,oerr=oerr,sigdig=sigdig,show.numbers=show.numbers,
                levels=levels,clabel=clabel,ellipse.fill=ellipse.fill,
                ellipse.stroke=ellipse.stroke,inverse=inverse,
                ci.col=ci.col,line.col=line.col,lwd=lwd,plot=plot,
                title=title,exterr=exterr,model=model,wtype=wtype,
                anchor=anchor,show.ellipses=show.ellipses,hide=hide,omit=omit,
                omit.fill=omit.fill,omit.stroke=omit.stroke,taxis=taxis,...)
}
#' @rdname isochron
#' @export
isochron.UThHe <- function(x,sigdig=2,oerr=3,show.numbers=FALSE,levels=NULL,
                           clabel="",ellipse.fill=c("#00FF0080","#FF000080"),
                           ellipse.stroke='black',ci.col='gray80',
                           line.col='black',lwd=1,plot=TRUE,title=TRUE,
                           model=1,wtype=1,anchor=0,show.ellipses=2*(model!=2),
                           hide=NULL,omit=NULL,omit.fill=NA,
                           omit.stroke='grey',...){
    wtype <- checkWtype(wtype=wtype,anchor=anchor,model=model)
    y <- data2york(x)
    d2calc <- clear(y,hide,omit)
    abanchor <- anchor
    if (anchor[1]==2 && length(anchor)>1){
        l2 <- lambda('Th232')[1]
        l5 <- lambda('U235')[1]
        l8 <- lambda('U238')[1]
        U <- iratio('U238U235')[1]
        He <- get_He(tt=anchor[2],U=1,Th=1)
        P <- 8*l8*U/(1+U) + 7*l5/(1+U) + 6*l2
        abanchor[2] <- He/P
        if (length(anchor)>2){
            abanchor[3] <- get_He(tt=anchor[3],U=1,Th=1)/P
        }
    }
    fit <- MLyork(d2calc,model=model,wtype=wtype,anchor=abanchor)
    out <- ab2y0t(x,fit=fit,wtype=wtype)
    if (plot){
        scatterplot(y,oerr=oerr,show.numbers=show.numbers,levels=levels,
                    clabel=clabel,ellipse.fill=ellipse.fill,
                    ellipse.stroke=ellipse.stroke,fit=out,
                    show.ellipses=show.ellipses,ci.col=ci.col,
                    line.col=line.col,lwd=lwd,hide=hide,omit=omit,
                    omit.fill=omit.fill,omit.stroke=omit.stroke,...)
        showDispersion(out,inverse=FALSE,wtype=wtype)
        dispunits <- getDispUnits(model=model,wtype=wtype,anchor=anchor)
        if (title){
            main <- isochrontitle(out,sigdig=sigdig,oerr=oerr,
                                  dispunits=dispunits,type='U-Th-He')
        } else {
            main <- NULL
        }
        graphics::title(main=main,xlab="P",ylab="He")
    }
    invisible(out)
}
#' @rdname isochron
#' @export
isochron.ThU <- function (x,type=2,oerr=3,sigdig=2,
                          show.numbers=FALSE,levels=NULL,clabel="",
                          ellipse.fill=c("#00FF0080","#FF000080"),
                          ellipse.stroke='black',ci.col='gray80',
                          line.col='black',lwd=1,plot=TRUE,
                          title=TRUE,exterr=FALSE,model=1,wtype='a',
                          show.ellipses=1*(model!=2),
                          hide=NULL,omit=NULL,omit.fill=NA,
                          omit.stroke='grey',y0option=4,...){
    displabel <- 'dispersion = '
    dispunits <- ''
    if (x$format %in% c(1,2)){
        out <- isochron_ThU_3D(x,type=type,model=model,wtype=wtype,
                               exterr=exterr,hide=hide,omit=omit,
                               y0option=y0option)
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
    } else if (x$format %in% c(3,4)){
        out <- isochron_ThU_2D(x,type=type,model=model,wtype=wtype,
                               exterr=exterr,hide=hide,omit=omit)
        if (model==3){
            if ((type==1 & wtype%in%c('slope',2,'b')) |
                (type==2 & wtype%in%c('intercept',1,'a'))){
                dispunits <- ' ka'
            } else {
                displabel <- quote('('^230*'Th/'^232*'Th)'[0]*'-dispersion = ')
            }
        }
    } else {
        stop('Illegal Th-U data format.')
    }
    out$disp <- w2disp(x,fit=out,type=type,wtype=wtype)
    if (plot){
        scatterplot(out$xyz,oerr=oerr,show.numbers=show.numbers,
                    levels=levels,clabel=clabel,ellipse.fill=ellipse.fill,
                    ellipse.stroke=ellipse.stroke,fit=out,
                    show.ellipses=show.ellipses,ci.col=ci.col,
                    line.col=line.col,lwd=lwd,hide=hide,omit=omit,
                    omit.fill=omit.fill,omit.stroke=omit.stroke,...)
        if (title){
            main <- isochrontitle(out,oerr=oerr,sigdig=sigdig,
                                  type='Th-U',units=' ka',
                                  displabel=displabel,dispunits=dispunits)
        } else {
            main <- NULL
        }
        graphics::title(main=main,xlab=out$xlab,ylab=out$ylab)
    }
    invisible(out)
}

isochron_ThU_2D <- function(x,type=2,model=1,wtype='a',
                            exterr=FALSE,hide=NULL,omit=NULL){
    yd <- data2york(x,type=type)
    d2calc <- clear(yd,hide,omit)
    out <- regression(d2calc,model=model,type="york",wtype=wtype)
    out$xyz <- yd
    if (type==1){
        Th230U238 <- out$b
        Th230Th232 <- out$a
    } else if (type==2) {
        Th230U238 <- out$a
        Th230Th232 <- out$b
    }
    cov0802 <- out$cov.ab
    
    l0 <- lambda('Th230')
    tt <- -log(1-Th230U238[1])/l0[1]
    y0 <- Th230Th232[1]/(1-Th230U238[1])
    J <- matrix(0,2,3)
    J[1,2] <- 1/(l0[1]*(1-Th230U238[1]))
    J[2,1] <- 1/(1-Th230U238[1])
    J[2,2] <- Th230Th232[1]/(1-Th230U238[1])^2
    if (exterr) J[1,3] <- -tt/l0[1]
    E <- matrix(0,3,3)
    diag(E) <- c(Th230Th232[2],Th230U238[2],l0[2])^2
    E[1,2] <- E[2,1] <- cov0802
    
    covmat = J %*% E %*% t(J)
    out$age[c('t','s[t]')] <- c(tt,sqrt(covmat[1,1]))
    out$y0[c('y','s[y]')] <- c(y0,sqrt(covmat[2,2]))
    
    if (inflate(out)){
        E[1:2,1:2] <- out$mswd*E[1:2,1:2]
        covmat = J %*% E %*% t(J)
        out$age['disp[t]'] <- sqrt(covmat[1,1])
        out$y0['disp[y]'] <- sqrt(covmat[2,2])
    }
    
    lab <- getIsochronLabels(x=x,type=type)
    out$xlab <- lab$x
    out$ylab <- lab$y
    out$y0label <- lab$y0
    out
}
isochron_ThU_3D <- function(x,type=2,model=1,wtype='a',exterr=FALSE,
                            hide=NULL,omit=NULL,y0option=4){
    if (type == 1){ # 0/2 vs 8/2
        osmond <- FALSE
        ia <- 'A'
        ib <- 'B'
        i48 <- 'b'
        i08 <- 'B'
        i02 <- 'A'
        id <- c('X','sX','Z','sZ','rXZ')
    } else if (type == 2){ # 0/8 vs 2/8
        osmond <- TRUE
        ia <- 'A'
        ib <- 'B'
        i48 <- 'a'
        i08 <- 'A'
        i02 <- 'B'
        id <- c('X','sX','Z','sZ','rXZ')
    } else if (type == 3){ # 4/2 vs 8/2
        osmond <- FALSE
        ia <- 'a'
        ib <- 'b'
        i48 <- 'b'
        i08 <- 'B'
        i02 <- 'A'
        id <- c('X','sX','Y','sY','rXY')
    } else if (type == 4){ # 4/8 vs 2/8
        osmond <- TRUE
        ia <- 'a'
        ib <- 'b'
        i48 <- 'a'
        i08 <- 'A'
        i02 <- 'B'
        id <- c('X','sX','Y','sY','rXY')
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
    tst <- get_ThU_age(out$par[i08],sqrt(out$cov[i08,i08]),
                       out$par[i48],sqrt(out$cov[i48,i48]),
                       out$cov[i48,i08],exterr=exterr,jacobian=TRUE)
    out$age['t'] <- tst['t']
    out$age['s[t]'] <- tst['s[t]']
    if (inflate(out)){
        tst <- get_ThU_age(out$par[i08],sqrt(out$mswd*out$cov[i08,i08]),
                           out$par[i48],sqrt(out$mswd*out$cov[i48,i48]),
                           out$mswd*out$cov[i48,i08],exterr=exterr,jacobian=TRUE)
        out$age['disp[t]'] <- tst['s[t]']
    }
    out <- getThUy0(out,tst=tst,option=y0option,exterr=exterr)
    lab <- getIsochronLabels(x=x,type=type)
    out$xlab <- lab$x
    out$ylab <- lab$y
    out$xyz <- subset(out$xyz,select=id)
    out
}

isochron_PD <- function(x,oerr=3,sigdig=2,show.numbers=FALSE,levels=NULL,
                        clabel="",ellipse.fill=c("#00FF0080","#FF000080"),
                        ellipse.stroke='black',inverse=FALSE,
                        ci.col='gray80',line.col='black',lwd=1,
                        plot=TRUE,title=TRUE,exterr=FALSE,model=1,
                        wtype=1,anchor=0,show.ellipses=1*(model!=2),
                        bratio=1,hide=NULL,omit=NULL,taxis=FALSE,...){
    taxis <- taxis & inverse
    wtype <- checkWtype(wtype=wtype,anchor=anchor,model=model)
    fit <- flipper(x,inverse=inverse,model=model,wtype=wtype,
                   anchor=anchor,hide=hide,omit=omit,type="p")
    out <- ab2y0t(x=x,fit=fit,inverse=inverse,wtype=wtype,
                  exterr=exterr,bratio=bratio)
    dispunits <- getDispUnits(model=model,wtype=wtype,anchor=anchor)
    lab <- getIsochronLabels(x=x,inverse=inverse,taxis=taxis)
    out$y0label <- lab$y0
    if (plot){
        scatterplot(out$xyz,oerr=oerr,show.ellipses=show.ellipses,
                    show.numbers=show.numbers,levels=levels,
                    clabel=clabel,ellipse.fill=ellipse.fill,
                    ellipse.stroke=ellipse.stroke,fit=out,
                    ci.col=ci.col,line.col=line.col,lwd=lwd,
                    hide=hide,omit=omit,taxis=taxis,...)
        if (taxis) add_taxis(x=x,fit=out,bratio=bratio)
        showDispersion(out,inverse=inverse,wtype=wtype)
        if (title){
            main <- isochrontitle(out,oerr=oerr,sigdig=sigdig,
                                  type='PD',dispunits=dispunits)
        } else {
            main <- NULL
        }
        graphics::title(main=main,xlab=lab$x,ylab=lab$y)
    }
    invisible(out)
}

get_isochron_PD_age <- function(DP,sDP,nuclide,exterr=FALSE,bratio=1,d=diseq()){
    if (nuclide=='U238'){
        out <- get_Pb206U238_age(x=DP,sx=sDP,exterr=exterr,d=d)
    } else if (nuclide=='U235'){
        out <- get_Pb207U235_age(x=DP,sx=sDP,exterr=exterr,d=d)
    } else {
        out <- getPDage(DP=DP,sDP=sDP,nuclide=nuclide,exterr=exterr,bratio=bratio)
    }
    out
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
    } else if (type=='generic'){
        content[[1]] <- maintit(x=fit$Dd[1],sx=fit$Dd[-1],n=fit$n,units=units,
                                prefix=quote('[D/d]'[0]*' ='),
                                sigdig=sigdig,oerr=oerr,df=fit$df)
        content[[2]] <- maintit(x=fit$DP[1],sx=fit$DP[-1],ntit='',units=units,
                                prefix=quote('[D/P]* ='),
                                sigdig=sigdig,oerr=oerr,df=fit$df)
    } else if (type=='U-Pb'){
        if (is.null(fit$posterior) | 't'%ni%names(fit$posterior)){
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
        else if (y0option==2 & 'U48i'%in%pnames) ipar <- 'U48i'
        else if (y0option==3 & 'ThUi'%in%pnames) ipar <-'ThUi'
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

getDispUnits_UPb <- function(x,joint,anchor){
    ifelse(anchor[1]==1 & (x$format%ni%(4:8) | !joint),'',' Ma')
}
getDispUnits <- function(model,wtype,anchor){
    ifelse(model==3 & (wtype==2 | anchor[1]==2), ' Ma','')
}

showDispersion <- function(fit,inverse,wtype,type='p'){
    if (fit$model!=3) return()
    usr <- graphics::par('usr')
    if (usr[1]>0 & usr[3]>0) return() # axes out of focus
    if (type=='p' & wtype==1){
        reldisp <- ifelse(inverse,
                          fit$flippedfit$w[1]/fit$flippedfit$a[1],
                          fit$w[1]/fit$a[1])
        y0 <- fit$a[1]
        cid <- ci(sx=y0*reldisp)
        graphics::lines(x=c(0,0),y=y0+cid*c(-1,1),lwd=2)
    } else if (type=='p' & wtype==2 & inverse){
        reldisp <- fit$flippedfit$w[1]/fit$flippedfit$a[1]
        x0 <- 1/fit$flippedfit$a[1]
        cid <- ci(sx=x0*reldisp)
        graphics::lines(x=x0+cid*c(-1,1),y=c(0,0),lwd=2)
    } else if (type=='d' & ((wtype==2 & inverse) | (wtype==1 & !inverse))){
        y0 <- fit$a[1]
        cid <- ci(sx=y0*fit$w[1]/fit$a[1])
        graphics::lines(x=c(0,0),y=y0+cid*c(-1,1),lwd=2)
    } else if (type=='TW' & wtype==1){
        y0 <- fit$par['a0']
        cid <- ci(sx=fit$disp[1])
        graphics::lines(x=c(0,0),y=y0+cid*c(-1,1),lwd=2)
    } else if (type%in%(1:4)){ # other UPb
        if (wtype==1){
            if (type==1){
                y0 <- 1/fit$par['a0']
                cid <- ci(sx=fit$disp[1]/fit$par['a0']^2)
            } else if (type==2){
                y0 <- 1/fit$par['b0']
                cid <- ci(sx=fit$disp[1]/fit$par['b0']^2)
            } else if (type==3){
                y0 <- fit$par['a0']
                cid <- ci(sx=fit$disp[1])
            } else if (type==4){
                y0 <- fit$par['b0']
                cid <- ci(sx=fit$disp[1])
            }
            graphics::lines(x=c(0,0),y=y0+cid*c(-1,1),lwd=2)
        } else {
            x0 <- -fit$a[1]/fit$b[1]
            cid <- ci(sx=x0*fit$disp[1]/fit$age[1])
            graphics::lines(x=x0+cid*c(-1,1),y=c(0,0),lwd=2)
        }
    }
}

#' Convert w parameter to meaningful dispersion estimate
#' @param x IsoplotR data object
#' @noRd
w2disp <- function(x,...){ UseMethod("w2disp",x) }
#' @param fit a linear model with a parameter w
#' @param wtype the dispersion type (0, 1, or 2)
#' @param inverse logical
#' @noRd
w2disp.default <- function(x,fit,wtype,inverse,...){ # type = 'p'
    if (wtype==2){
        out <- wDP2wt(x=x,DP=fit$flippedfit$a[1],
                      wDP=fit$flippedfit$w,...)
    } else if (inverse){
        out <- fit$y0[1]*fit$flippedfit$w/fit$flippedfit$a[1]
    } else {
        out <- fit$y0[1]*fit$w/fit$a[1]
    }
    out
}
#' @noRd
w2disp.other <- function(x,fit,wtype,inverse,...){
    if (wtype==2){
        out <- fit$DP[1]*fit$flippedfit$w/fit$flippedfit$a[1]
    } else if (inverse){
        out <- fit$Dd[1]*fit$flippedfit$w/fit$flippedfit$a[1]
    } else {
        out <- fit$Dd[1]*fit$w/fit$a[1]
    }
    out
}
#' @noRd
w2disp.PbPb <- function(x,fit,wtype,inverse,...){ # type = 'd'
    if (inverse){
        if (wtype==1){
            out <- fit$y0[1]*fit$flippedfit$w/fit$flippedfit$a[1]
        } else {
            out <- wDP2wt(x=x,DP=fit$a[1],wDP=fit$w)
        }
    } else {
        if (wtype==1){
            out <- fit$w
        } else {
            out <- wDP2wt(x=x,DP=fit$flippedfit$a[1],
                          wDP=fit$flippedfit$w)
        }
    }
    out
}
#' @noRd
w2disp.ThU <- function(x,fit,type,wtype,...){
    if (x$format%in%c(1,2)){
        out <- fit$w
    } else if (x$format%in%c(3,4)){
        if (type==1 & wtype%in%c('slope',2,'b')){
            age2disp <- TRUE
            Th230U238 <- fit$b[1]
        } else if (type==2 & wtype%in%c('intercept',1,'a')){
            age2disp <- TRUE
            Th230U238 <- fit$a[1]
        } else {
            age2disp <- FALSE
        }
        if (age2disp){
            out <- wDP2wt(x=x,DP=Th230U238,wDP=fit$w)
        } else {
            out <- fit$w
        }
    }
    out
}
#' @noRd
w2disp.UThHe <- function(x,fit,wtype,...){
    if (wtype==2){
        out <- fit$age[1]*fit$w/fit$b[1]
    } else {
        out <- fit$y0[1]*fit$w/fit$a[1]
    }
    out
}

#' Convert isotope ratio dispersion to age dispersion
#' @param x an IsoplotR data object
#' @noRd
wDP2wt <- function(x,...){ UseMethod("wDP2wt",x) }
#' @noRd
wDP2wt.default <- function(x,...){stop("Not implemented.")}
#' @noRd
wDP2wt.ArAr <- function(x,DP,wDP,...){
    l40 <- lambda('K40')[1]
    dtdDP <- x$J[1]/(l40*(1+DP*x$J[1]))
    abs(dtdDP*wDP)
}
#' @param DP a daughter-parent ratio
#' @param wDP the overdispersion of DP
#' @noRd
wDP2wt.ThPb <- function(x,DP,wDP,...){
    wDP2wt.PD(x=x,DP=DP,wDP=wDP,nuclide='Th232')
}
#' @param bratio branching ratio
#' @noRd
wDP2wt.KCa <- function(x,DP,wDP,bratio=0.895,...){
    wDP2wt.PD(x=x,DP=DP,wDP=wDP,nuclide='K40',bratio=bratio)
}
#' @param nuclide the parent nuclide
#' @noRd
wDP2wt.PD <- function(x,DP,wDP,nuclide,bratio=1,...){
    lambda <- lambda(nuclide)[1]
    dtdDP <- 1/(lambda*(1+DP*bratio))
    abs(dtdDP*wDP)
}
#' @noRd
wDP2wt.PbPb <- function(x,DP,wDP,...){
    tt <- get_Pb207Pb206_age(DP)[1]
    McL <- mclean(tt=tt)
    dtdDP <- 1/McL$dPb207Pb206dt
    abs(dtdDP*wDP)
}
#' @noRd
wDP2wt.ThU <- function(x,DP,wDP,...){
    l0 <- lambda('Th230')[1]
    dtdDP <- 1/(l0*(1-DP))
    abs(dtdDP*wDP)
}

#' Convert generic intercept and slope to inherited ratio and age
#' @param x an IsoplotR data object
#' @noRd
ab2y0t <- function(x,...){ UseMethod("ab2y0t",x) }
#' @param fit the output of york(), MLyork() or ludwig()
#' @noRd
ab2y0t.default <- function(x,fit,...){
    out <- fit
    if (inflate(fit)){
        out$a['disp[a]'] <- sqrt(fit$mswd)*fit$a['s[a]']
        out$b['disp[b]'] <- sqrt(fit$mswd)*fit$b['s[b]']
    } else if (fit$model==3){
        out$disp <- fit$w
    }
    out
}
#' @param type isochron type (see ?isochron)
#' @param exterr propagate decay constant errors?
#' @param y0option for datasets containing 204Pb or 208Pb (see ?isochron)
#' @noRd
ab2y0t.UPb <- function(x,fit,type,exterr=FALSE,y0option=1,...){
    tt <- fit$par['t']
    a0 <- fit$par['a0']
    b0 <- fit$par['b0']
    l8 <- lambda('U238')[1]
    l5 <- lambda('U235')[1]
    md <- mediand(x$d)
    if (md$U48$option==2) md$U48 <- list(x=unname(fit$par['U48i']),option=1)
    if (md$ThU$option==2) md$ThU <- list(x=unname(fit$par['ThUi']),option=1)
    McL <- mclean(tt,d=md,exterr=exterr)
    if (type==1){                           # 04-08c/06 vs. 38/06
        x0inv <- McL$Pb206U238
        dx0invdt <- McL$dPb206U238dt
        E <- fit$cov[c('t','a0'),c('t','a0')]
    } else if (type==2){                    # 04-08c/07 vs. 35/07
        x0inv <- McL$Pb207U235
        dx0invdt <- McL$dPb207U235dt
        E <- fit$cov[c('t','b0'),c('t','b0')]
    } else if (type==3 & x$format%in%c(7,8)){  # 06c/08 vs. 32/08
        x0inv <- age_to_Pb208Th232_ratio(tt=tt,st=0)[1]
        dx0invdt <- McL$dPb208Th232dt
        E <- fit$cov[c('t','a0'),c('t','a0')]
    } else if (type==4 & x$format%in%c(7,8)){  # 07c/08 vs. 32/08
        x0inv <- age_to_Pb208Th232_ratio(tt=tt,st=0)[1]
        dx0invdt <- McL$dPb208Th232dt
        E <- fit$cov[c('t','b0'),c('t','b0')]
    } else {
        stop('Invalid isochron type.')
    }
    out <- fit
    out$age <- NULL
    out$age['t'] <- tt
    out$age['s[t]'] <- sqrt(fit$cov['t','t'])
    J <- matrix(0,2,2)
    if (x$format%in%c(4,5,6,9,85,119) & type==1){          # 0x/06 vs. 38/06
        out$XY <- data2york(x,option=3)
        a <- 1/fit$par['a0']
        J[1,2] <- -a^2
    } else if (x$format%in%c(4,5,6,10,85,1210) & type==2){ # 0x/07 vs. 35/07
        out$XY <- data2york(x,option=4)
        a <- 1/fit$par['b0']
        J[1,2] <- -a^2
    } else if (x$format%in%c(7,8,11) & type==1){   # 08c/06 vs. 38/06
        out$XY <- data2york(x,option=6,tt=tt)
        a <- 1/fit$par['a0']
        J[1,2] <- -a^2
    } else if (x$format%in%c(7,8,12) & type==2){   # 08c/07 vs. 35/07
        out$XY <- data2york(x,option=7,tt=tt)
        U <- settings('iratio','U238U235')[1]
        a <- 1/fit$par['b0']
        J[1,2] <- -a^2
    } else if (x$format%in%c(7,8,11) & type==3){   # 06c/08 vs. 32/08
        out$XY <- data2york(x,option=8,tt=tt)
        a <- fit$par['a0']
        J[1,2] <- 1
    } else if (x$format%in%c(7,8,12) & type==4){   # 07c/08 vs. 32/08
        out$XY <- data2york(x,option=9,tt=tt)
        a <- fit$par['b0']
        J[1,2] <- 1
    } else {
        stop('Isochron regression is not available for this input format.')
    }
    b <- -a*x0inv
    J[2,1] <- -a*dx0invdt
    J[2,2] <- x0inv*a^2
    cov.ab <- J%*%E%*%t(J)
    out$a <- c(a,sqrt(cov.ab[1,1]))
    out$b <- c(b,sqrt(cov.ab[2,2]))
    names(out$a) <- c('a','s[a]')
    names(out$b) <- c('b','s[b]')
    out$cov.ab <- cov.ab[1,2]
    out <- getUPby0(out,fmt=x$format,type=type,option=y0option)
    if (inflate(out)){
        out$age['disp[t]'] <- sqrt(out$mswd)*out$age['s[t]']
        out$y0['disp[y]'] <- sqrt(out$mswd)*out$y0['s[y]']
    } else if (out$model==3){ # wx0 and wy0 set error bars in scatterplot
        out$wx0 <- out$wy0 <- 0
        if (fit$wtype==1){
            out$wy0 <- out$disp[1]
        } else if (type==1){ # 06/38
            Pb6U8 <- age2ratio(out$age[1],out$disp[1],ratio='Pb206U238',d=x$d)
            out$wx0 <- Pb6U8[2]/Pb6U8[1]^2
        } else if (type==2){ # 07/35
            Pb7U5 <- age2ratio(out$age[1],out$disp[1],ratio='Pb207U235',d=x$d)
            out$wx0 <- Pb7U5[2]/Pb7U5[1]^2
        } else if (type%in%c(3,4)){ # 08/32
            Pb8Th2 <- age2ratio(out$age[1],out$disp[1],ratio='Pb208Th232')
            out$wx0 <- Pb8Th2[2]/Pb8Th2[1]^2
        } else {
            stop('Invalid isochron type')
        }
    }
    out
}
#' @param inverse logical
#' @param wtype controls type of dispersion (0, 1 or 2)
#' @noRd
ab2y0t.PbPb <- function(x,fit,inverse,exterr,wtype,...){
    out <- fit
    if (inverse){
        R76 <- fit$a
        out$y0[c('y','s[y]')] <- fit$b
    } else {
        R76 <- fit$b
        out$y0[c('y','s[y]')] <- fit$a
    }
    out$y0label <- quote('('^207*'Pb/'^204*'Pb)'[c]*'=')
    out$age[c('t','s[t]')] <- get_Pb207Pb206_age(R76[1],R76[2],exterr=exterr)
    if (inflate(fit)){
        out$age['disp[t]'] <- 
            get_Pb207Pb206_age(R76[1],sqrt(fit$mswd)*R76[2],exterr=exterr)[2]
        out$y0['disp[y]'] <- sqrt(fit$mswd)*out$y0['s[y]']
    } else if (fit$model==3){
        out$disp <- w2disp(x=x,fit=out,wtype=wtype,inverse=inverse)
    }
    out
}
#' @noRd
ab2y0t.ArAr <- function(x,fit,inverse,exterr,wtype,...){
    out <- fit
    if (inverse) {
        R09 <- quotient(X=fit$a[1],sX=fit$a[2],
                        Y=fit$b[1],sY=fit$b[2],sXY=fit$cov.ab)
        R09[1] <- -R09[1]
        out$y0[c('y','s[y]')] <- quotient(X=fit$a[1],sX=fit$a[2],Y=1,sY=0,rXY=0)
    } else {
        R09 <- fit$b
        out$y0[c('y','s[y]')] <- fit$a
    }
    out$y0label <- quote('('^40*'Ar/'^36*'Ar)'[0]*'=')
    out$age[c('t','s[t]')] <- get_ArAr_age(R09[1],R09[2],x$J[1],x$J[2],exterr=exterr)
    if (inflate(out)){
        out$age['disp[t]'] <- get_ArAr_age(R09[1],sqrt(fit$mswd)*R09[2],
                                           x$J[1],x$J[2],exterr=exterr)[2]
        out$y0['disp[y]'] <- sqrt(fit$mswd)*out$y0['s[y]']
    } else if (fit$model==3){
        out$disp <- w2disp(x=x,fit=out,wtype=wtype,inverse=inverse)
    }
    out
}
#' @noRd
ab2y0t.ThPb <- function(x,fit,inverse,exterr,wtype,...){
    ab2y0t.PD(x=x,fit=fit,inverse=inverse,exterr=exterr,
              nuclide='Th232',wtype=wtype)
}
#' @noRd
ab2y0t.KCa <- function(x,fit,inverse,exterr,bratio=0.895,wtype,...){
    ab2y0t.PD(x=x,fit=fit,inverse=inverse,exterr=exterr,
              nuclide='K40',bratio=bratio,wtype=wtype)
}
#' @param bratio branching ratio
#' @noRd
ab2y0t.PD <- function(x,fit,inverse,exterr,bratio=1,wtype,...){
    nuclide <- getParent(x)
    out <- fit
    if (inverse){
        Dd <- c(1,fit$a[2]/fit$a[1])/fit$a[1]
        DP <- quotient(X=fit$a[1],sX=fit$a[2],
                       Y=fit$b[1],sY=fit$b[2],sXY=fit$cov.ab)
        DP[1] <- -DP[1]
    } else {
        Dd <- fit$a
        DP <- fit$b
    }
    out$y0[c('y','s[y]')] <- unname(Dd)
    out$age[c('t','s[t]')] <-
        get_isochron_PD_age(DP=DP[1],sDP=DP[2],nuclide=nuclide,
                            exterr=exterr,bratio=bratio)
    if (inflate(fit)){
        out$age['disp[t]'] <-
            get_isochron_PD_age(DP=DP[1],sDP=sqrt(fit$mswd)*DP[2],
                                nuclide=nuclide,exterr=exterr,bratio=bratio)[2]
        out$y0['disp[y]'] <- sqrt(fit$mswd)*out$y0['s[y]']
    } else if (fit$model==3){
        out$disp <- w2disp(x=x,fit=out,wtype=wtype,inverse=inverse,nuclide=nuclide)
    }
    out
}
#' @noRd
ab2y0t.UThHe <- function(x,fit,wtype,...){
    out <- fit
    out$y0[c('y','s[y]')] <- fit$a
    out$age[c('t','s[t]')] <- fit$b
    if (inflate(out)){
        out$age['disp[t]'] <- sqrt(fit$mswd)*out$age['s[t]']
        out$y0['disp[y]'] <- sqrt(fit$mswd)*out$y0['s[y]']
    } else if (fit$model==3){
        out$disp <- w2disp(x=x,fit=out,wtype=wtype,inverse=FALSE)
    }
    out$y0label <- quote('He'[0]*'=')
    out
}
#' @noRd
ab2y0t.other <- function(x,fit,inverse,wtype,...){
    out <- fit
    if (inverse){
        Dd <- 1/fit$a[1]
        DP <- -fit$b[1]/fit$a[1]
        err <- errorprop(J11=-Dd/fit$a[1],
                         J12=0,
                         J21=-DP/fit$a[1],
                         J22=-1/fit$a[1],
                         E11=fit$a[2]^2,
                         E22=fit$b[2]^2,
                         E12=fit$cov.ab)
        out$Dd <- unname(c(Dd,sqrt(err[1])))
        out$DP <- unname(c(DP,sqrt(err[2])))
        out$cov.DdDp <- unname(err[3])
    } else {
        out$Dd <- unname(fit$a)
        out$DP <- unname(fit$b)
        out$cov.DdDP <- unname(fit$cov.ab)
    }
    names(out$Dd) <- c('Dd','s[Dd]')
    names(out$DP) <- c('DP','s[DP]')
    if (inflate(fit)){
        out$Dd['disp[Dd]'] <- sqrt(fit$mswd)*out$Dd['s[Dd]']
        out$DP['disp[DP]'] <- sqrt(fit$mswd)*out$DP['s[DP]']
    } else if (fit$model==3){
        out$disp <- w2disp(x=x,fit=out,wtype=wtype,inverse=inverse)
    }
    out
}

#' Generates expressions for x- and y- axis labels of isochron plots
#' @param x an IsoplotR data object
#' @noRd
getIsochronLabels <- function(x,...){ UseMethod("getIsochronLabels",x) }
#' @noRd
getIsochronLabels.default <- function(x,...){stop("Not implemented.")}
#' @param inverse logical
#' @param taxis label the x-axis with time
#' @noRd
getIsochronLabels.ThPb <- function(x,inverse,taxis=FALSE,...){
    out <- list()
    if (taxis){
        out$x <- 't (Ma)'
    } else if (inverse){
        out$x <- quote(''^232*'Th/'^208*'Pb')
    } else {
        out$x <- quote(''^232*'Th/'^204*'Pb')
    }
    if (inverse){
        out$y <- quote(''^204*'Pb/'^208*'Pb')
    } else {
        out$y <- quote(''^208*'Pb/'^204*'Pb')
    }
    out$y0 <- quote('('^208*'Pb/'^204*'Pb)'[0]*'=')
    out
}
#' @noRd
getIsochronLabels.SmNd <- function(x,inverse,taxis=FALSE,...){
    out <- list()
    if (taxis){
        out$x <- 't (Ma)'
    } else if (inverse){
        out$x <- quote(''^147*'Sm/'^143*'Nd')
    } else {
        out$x <- quote(''^147*'Sm/'^144*'Nd')
    }
    if (inverse){
        out$y <- quote(''^144*'Nd/'^143*'Nd')
    } else {
        out$y <- quote(''^143*'Nd/'^144*'Nd')
    }
    out$y0 <- quote('('^143*'Nd/'^144*'Nd)'[0]*'=')
    out
}
#' @noRd
getIsochronLabels.ReOs <- function(x,inverse,taxis=FALSE,...){
    out <- list()
    if (taxis){
        out$x <- 't (Ma)'
    } else if (inverse){
        out$x <- quote(''^187*'Re/'^187*'Os')
    } else {
        out$x <- quote(''^187*'Re/'^188*'Os')
    }
    if (inverse){
        out$y <- quote(''^188*'Os/'^187*'Os')
    } else {
        out$y <- quote(''^187*'Os/'^188*'Os')
    }
    out$y0 <- quote('('^187*'Os/'^188*'Os)'[0]*'=')
    out
}
#' @noRd
getIsochronLabels.RbSr <- function(x,inverse,taxis=FALSE,...){
    out <- list()
    if (taxis){
        out$x <- 't (Ma)'
    } else if (inverse){
        out$x <- quote(''^87*'Rb/'^87*'Sr')
    } else {
        out$x <- quote(''^87*'Rb/'^86*'Sr')
    }
    if (inverse){
        out$y <- quote(''^86*'Sr/'^87*'Sr')
    } else {
        out$y <- quote(''^87*'Sr/'^86*'Sr')
    }
    out$y0 <- quote('('^87*'Sr/'^86*'Sr)'[0]*'=')
    out
}
#' @noRd
getIsochronLabels.LuHf <- function(x,inverse,taxis=FALSE,...){
    out <- list()
    if (taxis){
        out$x <- 't (Ma)'
    } else if (inverse){
        out$x <- quote(''^176*'Lu/'^176*'Hf')
    } else {
        out$x <- quote(''^176*'Lu/'^177*'Hf')
    }
    if (inverse){
        out$y <- quote(''^177*'Hf/'^176*'Hf')
    } else {
        out$y <- quote(''^176*'Hf/'^177*'Hf')
    }
    out$y0 <- quote('('^176*'Hf/'^177*'Hf)'[0]*'=')
    out
}
#' @noRd
getIsochronLabels.KCa <- function(x,inverse,taxis=FALSE,...){
    out <- list()
    if (taxis){
        out$x <- 't (Ma)'
    } else if (inverse){
        out$x <- quote(''^40*'K/'^40*'Ca')
    } else {
        out$x <- substitute(''^40*'K/'^s*'Ca',list(s=x$sister))
    }
    if (inverse){
        out$y <- substitute(''^s*'Ca/'^40*'Ca',list(s=x$sister))
    } else {
        out$y <- substitute(''^40*'Ca/'^s*'Ca',list(s=x$sister))
    }
    out$y0 <- substitute('('^40*'Ca/'^s*'Ca)'[0]*'=',list(s=x$sister))
    out
}
#' @param type controls the isochron projection
#' @noRd
getIsochronLabels.UPb <- function(x,type=1,taxis=FALSE,...){
    out <- list()
    if (taxis){
        out$x <- 't (Ma)'
    } else if (type==1){
        out$x <- quote(''^238*'U/'^206*'Pb')
    } else if (type==2){
        out$x <- quote(''^235*'U/'^207*'Pb')
    } else if (type==3 & x$format%in%c(7,8)){
        out$x <- quote(''^232*'Th/'^208*'Pb')
    } else if (type==4 & x$format%in%c(7,8)){
        out$x <- quote(''^232*'Th/'^208*'Pb')
    } else {
        stop('Invalid isochron type.')
    }
    if (x$format%in%c(4,5,6,9) & type==1){
        out$y <- quote(''^204*'Pb/'^206*'Pb')
    } else if (x$format%in%c(85,119) & type==1){
        out$y <- quote(''^208*'Pb/'^206*'Pb')
    } else if (x$format%in%c(4,5,6,10) & type==2){
        out$y <- quote(''^204*'Pb/'^207*'Pb')
    } else if (x$format%in%c(85,1210) & type==2){
        out$y <- quote(''^208*'Pb/'^207*'Pb')
    } else if (x$format%in%c(7,8,11) & type==1){
        out$y <- quote(''^208*'Pb'[c]*'/'^206*'Pb')
    } else if (x$format%in%c(7,8,12) & type==2){
        out$y <- quote(''^208*'Pb'[c]*'/'^207*'Pb')
    } else if (x$format%in%c(7,8,11) & type==3){
        out$y <- quote(''^206*'Pb'[c]*'/'^208*'Pb')
    } else if (x$format%in%c(7,8,12) & type==4){
        out$y <- quote(''^207*'Pb'[c]*'/'^208*'Pb')
    } else {
        stop('Invalid isochron type.')
    }
    out
}
#' @noRd
getIsochronLabels.PbPb <- function(x,inverse,...){
    out <- list()
    if (inverse){
        out$x <- quote(''^204*'Pb/'^206*'Pb')
        out$y <- quote(''^207*'Pb/'^206*'Pb')
    } else {
        out$x <- quote(''^206*'Pb/'^204*'Pb')
        out$y <- quote(''^207*'Pb/'^204*'Pb')
    }
    out$y0 <- quote('('^207*'Pb/'^204*'Pb)'[c]*'=')
    out
}
#' @noRd
getIsochronLabels.ArAr <- function(x,inverse,taxis=FALSE,...){
    out <- list()
    if (taxis){
        out$x <- 't (Ma)'
    } else if (inverse){
        out$x <- quote(''^39*'Ar/'^40*'Ar')
    } else {
        out$x <- quote(''^39*'Ar/'^36*'Ar')
    }
    if (inverse){
        out$y <- quote(''^36*'Ar/'^40*'Ar')
    } else {
        out$y <- quote(''^40*'Ar/'^36*'Ar')
    }
    out$y0 <- quote('('^40*'Ar/'^36*'Ar)'[0]*'=')
    out
}
#' @noRd
getIsochronLabels.other <- function(x,inverse,...){
    out <- list()
    if (inverse){
        xlab <- 'P/D'
        ylab <- 'd/D'
    } else {
        xlab <- 'P/d'
        ylab <- 'D/d'
    }
    out$y0 <- quote('[D/d]'[0]*'=')
    out
}
#' @noRd
getIsochronLabels.ThU <- function(x,type,...){
    out <- list()
    if (x$format%in%c(1,2)){
        if (type == 1){ # 0/2 vs 8/2
            out$x <- quote(''^238*'U/'^232*'Th')
            out$y <- quote(''^230*'Th/'^232*'Th')
        } else if (type == 2){ # 0/8 vs 2/8
            out$x <- quote(''^232*'Th/'^238*'U')
            out$y <- quote(''^230*'Th/'^238*'U')
        } else if (type == 3){ # 4/2 vs 8/2
            out$x <- quote(''^238*'U/'^232*'Th')
            out$y <- quote(''^234*'U/'^232*'Th')
        } else if (type == 4){ # 4/8 vs 2/8
            out$x <- quote(''^232*'Th/'^238*'U')
            out$y <- quote(''^234*'U/'^238*'U')
        }
    } else if (x$format%in%c(3,4)){
        if (type==1){
            out$x <- quote(''^238*'U/'^232*'Th')
            out$y <- quote(''^230*'Th/'^232*'Th')
        } else if (type==2) {
            out$x <- quote(''^232*'Th/'^238*'U')
            out$y <- quote(''^230*'Th/'^238*'U')
        }
        out$y0label <- quote('('^230*'Th/'^232*'Th)'[0]*'=')
    } else {
        stop('Invalid input format')
    }
    out
}

getUPby0 <- function(out,fmt=1,type=1,option=1){
    out$y0 <- c('y'=NA,'s[y]'=NA)
    if (option==1){
        if (fmt<4){                                       # 07/06 vs. 38/06
            out$y0['y'] <- out$par['a0']
            out$y0['s[y]'] <- out$err['s','a0']
            out$y0label <- quote('('^207*'Pb/'^206*'Pb)'[c]*'=')
        } else if (fmt %in% c(4,5,6,9,85,119) & type==1){ # 0x/06 vs. 38/06
            out$y0['y'] <- out$par['a0']
            out$y0['s[y]'] <- sqrt(out$cov['a0','a0'])
            if (fmt%in%c(85,119)){
                out$y0label <- quote('('^206*'Pb/'^208*'Pb)'[c]*'=')
            } else {
                out$y0label <- quote('('^206*'Pb/'^204*'Pb)'[c]*'=')
            }
        } else if (fmt %in% c(4,5,6,10,85,1210) & type==2){ # 0x/07 vs. 35/07
            out$y0['y'] <- out$par['b0']
            out$y0['s[y]'] <- sqrt(out$cov['b0','b0'])
            if (fmt%in%c(85,1210)){
                out$y0label <- quote('('^207*'Pb/'^208*'Pb)'[c]*'=')
            } else {
                out$y0label <- quote('('^207*'Pb/'^204*'Pb)'[c]*'=')
            }
        } else if (fmt %in% c(7,8,11) & type%in%c(1,3)){
            out$y0['y'] <- out$par['a0']
            out$y0['s[y]'] <- sqrt(out$cov['a0','a0'])
            out$y0label <- quote('('^206*'Pb/'^208*'Pb)'[c]*'=')
        } else if (fmt %in% c(7,8,12) & type%in%c(2,4)){   # 08/07 vs. 35/07
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

#' Add a time axis to a scatterplot
#' @param x an IsoplotR data object
#' @noRd
add_taxis <- function(x,...){ UseMethod("add_taxis",x) }
#' @param fit the output of york(), MLyork() or ludwig()
#' @noRd
add_taxis.default <- function(x,fit,...){
    xlim <- graphics::par('usr')[1:2]
    xmid <- xlim[1] + diff(xlim)/3
    ratio <- getDPrat(x)
    parent <- getParent(x)
    tmin <- getPDage(DP=1/xlim[2],nuclide=parent,...)[1]
    xzero <- 1/age2ratio(tt=5000,st=0,ratio=ratio,...)[1]
    if (xzero<xmid){ # 5Ga is to the left of the middle
        tmid <- getPDage(DP=1/xmid,sDP=0,nuclide=parent,...)[1]
    } else {
        tmid <- getPDage(DP=1/xzero,sDP=0,nuclide=parent,...)[1]
    }
    plot_taxis(x=x,fit=fit,tmin=tmin,tmid=tmid,ratio=ratio,...)
}
#' @noRd
add_taxis.ArAr <- function(x,fit,...){
    xlim <- graphics::par('usr')[1:2]
    xmid <- xlim[1] + diff(xlim)/3
    ratio <- 'Ar40Ar39'
    tmin <- get_ArAr_age(Ar40Ar39=1/xlim[2],J=x$J[1])[1]
    xzero <- 1/age2ratio(tt=5000,ratio=ratio,J=x$J[1])[1]
    if (xzero<xmid){ # 5Ga is to the left of the middle
        tmid <- get_ArAr_age(Ar40Ar39=1/xmid,J=x$J[1])[1]
    } else {
        tmid <- get_ArAr_age(Ar40Ar39=1/xzero,J=x$J[1])[1]
    }
    plot_taxis(x=x,fit=fit,tmin=tmin,tmid=tmid,ratio=ratio,J=x$J[1])
}
#' @param type controls the isochron projection
#' @noRd
add_taxis.UPb <- function(x,fit,type=1,...){
    xlim <- graphics::par('usr')[1:2]
    xmid <- xlim[1] + diff(xlim)/3
    if (type==1){
        tfun <- get_Pb206U238_age
        ratio <- 'Pb206U238'
    } else if (type==2){
        tfun <- get_Pb207U235_age
        ratio <- 'Pb207U235'
    } else if (type==3){
        tfun <- get_Pb208Th232_age
        ratio <- 'Pb208Th232'
    } else {
        stop("Invalid type")
    }
    tmin <- tfun(x=1/xlim[2],d=x$d)[1]
    xzero <- 1/age2ratio(tt=5000,ratio=ratio,d=x$d)[1]
    if (xzero<xmid){ # 5Ga is to the left of the middle
        tmid <- tfun(x=1/xmid,sx=0,d=x$d)[1]
    } else {
        tmid <- tfun(x=1/xzero,sx=0,d=x$d)[1]
    }
    plot_taxis(x=x,fit=fit,tmin=tmin,tmid=tmid,ratio=ratio,d=x$d,...)
}

plot_taxis <- function(x,fit,tmin,tmid,ratio,...){
    tmax <- max(tmid,fit$age[1] + 3*utils::tail(fit$age,1))
    tlim <- c(tmin,tmax)
    tticks <- taxisticks(tlim)
    xticks <- 1/age2ratio(tt=tticks,ratio=ratio,...)[,1]
    graphics::axis(1,at=xticks,labels=signif(tticks,5))
}

taxisticks <- function(tlim){
    tticks <- pretty(tlim)
    extraticks <- pretty(tticks[1:2])
    tticks[1] <- min(extraticks[extraticks>min(tlim)])
    tticks
}
