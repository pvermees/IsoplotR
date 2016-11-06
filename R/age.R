#' Calculate isotopic ages
#'
#' Calculates  ages and propagates their analytical
#' uncertainties.
#'
#' @param x can be:
#'
#' - a scalar containing an isotopic ratio,
#'
#' - a two element vector containing an isotopic ratio and its standard
#' error, or the spontaneous and induced track densities \code{Ns} and
#' \code{Ni} (if \code{method='fissiontracks'}),
#'
#' - a four element vector containing \code{Ar40Ar39},
#' \code{s[Ar40Ar39]}, \code{J}, \code{s[J]},
#'
#' - a six element vector containing \code{U}, \code{s[U]}, \code{Th},
#' \code{s[Th]}, \code{He} and \code{s[He]},
#'
#' - an eight element vector containing \code{U}, \code{s[U]},
#' \code{Th}, \code{s[Th]}, \code{He}, \code{s[He]}, \code{Sm} and
#' \code{s[Sm]},
#' 
#' OR
#' 
#' - an object of class \code{UPb}, \code{ArAr}, \code{UThHe} or
#' \code{fissiontracks}.
#' 
#' @param method one of either \code{'Pb206U238'}, \code{'Pb207U235'},
#'     \code{'Pb207Pb206'}, \code{'Ar40Ar39'}, \code{U-Th-He} or
#'     \code{fissiontracks}
#' 
#' @param exterr propagate the external (decay constant and
#'     calibration factor) uncertainties?
#' 
#' @param J two-element vector with the J-factor and its standard
#'     error.  This option is only used if \code{method} =
#'     \code{'Ar40Ar39'}.
#' 
#' @param zeta two-element vector with the zeta-factor and its standard
#'     error.  This option is only used if \code{method} =
#'     \code{'fissiontracks'}.
#' 
#' @param rhoD two-element vector with the track density of the
#'     dosimeter glass and its standard error.  This option is only
#'     used if \code{method} = \code{'fissiontracks'}.
#' 
#' @param i (optional) index of a particular aliquot
#' 
#' @param ... optional arguments
#' 
#' @rdname age
#' @export
age <- function(x,...){ UseMethod("age",x) }
#' @rdname age
#' @export
age.default <- function(x,method='Pb206U238',exterr=TRUE,J=c(NA,NA),
                        zeta=c(NA,NA),rhoD=c(NA,NA),...){
    if (length(x)==1) X <- c(x,0)
    else X <- x[1:2]
    if (identical(method,'Pb207U235')){
        out <- get.Pb207U235age(X[1],X[2],exterr)
    } else if (identical(method,'Pb206U238')){
        out <- get.Pb206U238age(X[1],X[2],exterr)
    } else if (identical(method,'Pb207Pb206')){
        out <- get.Pb207Pb206age(X[1],X[2],exterr)
    } else if (identical(method,'Ar40Ar39')){
        out <- get.ArAr.age(X[1],X[2],X[3],X[4],exterr)
    } else if (identical(method,'U-Th-He')){
        if (length(x)==6)
            out <- get.UThHe.age(X[1],X[2],X[3],X[4],X[5],X[6])
        else if (length(x)==8)
            out <- get.UThHe.age(X[1],X[2],X[3],X[4],X[5],X[6],X[7],X[8])
    } else if (identical(method,'fissiontracks')){
        out <- get.EDM.age(X[1],X[2],zeta,rhoD)
    }
    out
}
#' @param concordia scalar flag indicating whether each U-Pb analysis
#'     should be considered separately (\code{concordia=1}), a
#'     concordia age should be calculated from all U-Pb analyses
#'     together (\code{concordia=2}), or a discordia line should be
#'     fit through all the U-Pb analyses (\code{concordia=3}).
#' 
#' @param wetherill logical flag to indicate whether the data should
#'     be evaluated in Wetherill (\code{TRUE}) or Tera-Wasserburg
#'     (\code{FALSE}) space.  This option is only used when
#'     \code{concordia=2}
#' 
#' @param sigdig number of significant digits for the uncertainty
#'     estimate (only used if \code{concordia=1}, \code{isochron=FALSE}
#'     or \code{central=FALSE}).
#'
#' @return 
#'
#' \enumerate{
#'
#' \item if \code{x} is a scalar or a vector, returns the age using
#' the geochronometer given by \code{method} and its standard error.
#'
#' \item if \code{x} has class \code{UPb} and \code{concordia=1},
#' returns a table with the following columns: `t.75', `err[t.75]',
#' `t.68', `err[t.68]', `t.76',`err[t.76]', `t.conc', `err[t.conc]',
#' containing the \eqn{^{207}}Pb/\eqn{^{235}} U-age and standard
#' error, the \eqn{^{206}}Pb/\eqn{^{238}}U-age and standard error, the
#' \eqn{^{207}}Pb/\eqn{^{206}} Pb-age and standard error, and the
#' concordia age and standard error, respectively.
#'  
#' \item if \code{x} has class \code{UPb} and \code{concordia=2},
#' returns a list with the following items:
#'
#' \describe{
#' \item{x}{ a named vector with the (weighted mean) U-Pb composition }
#' 
#' \item{cov}{ the covariance matrix of the (mean) U-Pb composition }
#' 
#' \item{age}{ the concordia age (in Ma) }
#' 
#' \item{age.err}{ the standard error of the concordia age }
#' 
#' \item{mswd}{ a list with two items (\code{equivalence} and
#' \code{concordance}) containing the MSWD (Mean of the Squared
#' Weighted Deviates, a.k.a the reduced Chi-squared statistic outside
#' of geochronology) of isotopic equivalence and age concordance,
#' respectively. }
#' 
#' \item{p.value}{ a list with two items (\code{equivalence} and
#' \code{concordance}) containing the p-value of the Chi-square test
#' for isotopic equivalence and age concordance, respectively. }
#' }
#' 
#' \item if \code{x} has class \code{UPb} and \code{concordia=3},
#' returns a list with the following items:
#'
#' \describe{
#' \item{x}{ a two element vector with the upper and lower intercept
#' ages (if wetherill==TRUE) or the lower intercept age and
#' \eqn{^{207}}Pb/\eqn{^{206}}Pb intercept (for Tera-Wasserburg) }
#' 
#' \item{cov}{ the covariance matrix of the elements in \code{x} }
#' }
#'
#' \item if \code{x} has class \code{ArAr} and \code{isochron=FALSE},
#' returns a table of Ar-Ar ages and standard errors.
#'
#' \item if \code{x} has class \code{ArAr} and \code{isochron=TRUE},
#' returns a list with the following items:
#'
#' \describe{
#'
#' \item{a}{ the intercept of the straight line fit and its standard
#' error. }
#' 
#' \item{b}{ the slope of the fit and its standard error. }
#' 
#' \item{y0}{ the atmospheric \eqn{^{40}}Ar/\eqn{^{36}}Ar ratio and
#' its standard error. }
#' 
#' \item{age}{ the \eqn{^{40}}Ar/\eqn{^{39}}Ar age and its standard
#' error. }
#' 
#' }
#' 
#' \item if \code{x} has class \code{UThHe} and \code{central=FALSE},
#' returns a table of U-Th-He ages and standard errors.
#' 
#' \item if \code{x} has class \code{UThHe} and \code{central=TRUE},
#' returns a list with the following items:
#'
#' \describe{
#'
#' \item{uvw}{ a three-element list with the weighted mean log[U/He],
#' log[Th/He] and log[Sm/He] compositions. }
#'
#' \item{covmat}{ a 3x3 covariance matrix for uvw}
#'
#' \item{mswd}{ the reduced Chi-square value for the
#' log[U/He]-log[Th/He] compositions. }
#'
#' \item{p.value}{ the p-value of concordance between the
#' log[U/He]-log[Th/He] compositions. }
#'
#' \item{age}{ two-element vector with the central age and its
#' standard error. }
#'
#' }
#'
#' \item if \code{x} has class \code{fissiontracks} and
#' \code{central=FALSE}, returns a table of fission track ages and
#' standard errors.
#' 
#' \item if \code{x} has class \code{fissiontracks} and
#' \code{central=TRUE}, returns a list with the following items:
#'
#' \describe{
#'
#' \item{mswd}{ the reduced Chi-square value for the fission track
#' ages. }
#'
#' \item{p.value}{ the p-value of concordance between the fission
#' track ages. }
#'
#' \item{age}{ a two-element vector with the central age and its
#' standard error. }
#'
#' \item{disp}{ the (over)dispersion of the single grain ages beyond
#' the formal analytical uncertainties. }
#' }
#' }
#' 
#' @examples
#' data(examples)
#' print(age(examples$UPb))
#' print(age(examples$UPb,concordia=1))
#' print(age(examples$UPb,concordia=2))
#' @rdname age
#' @export
age.UPb <- function(x,concordia=1,wetherill=TRUE,
                    exterr=TRUE,i=NA,sigdig=NA,...){
    if (concordia==1)
        out <- UPb.age(x,exterr=exterr,i=i,sigdig=sigdig,...)
    else if (concordia==2)
        out <- concordia.age(x,wetherill=TRUE,exterr=TRUE,...)
    else if (concordia==3)
        out <- discordia.age(x,wetherill=TRUE,exterr=TRUE,...)
    out
}
#' @rdname age
#' @export
age.detritals <- function(x,...){
    x
}
#' @param isochron logical flag indicating whether each Ar-Ar analysis
#'     should be considered separately (\code{isochron=FALSE}) or an
#'     isochron age should be calculated from all Ar-Ar analyses
#'     together (\code{isochron=TRUE}).
#' @rdname age
#' @export
age.ArAr <- function(x,isochron=FALSE,exterr=TRUE,i=NA,sigdig=NA,...){
    if (isochron) out <- isochron(x,plot=FALSE)
    else out <- ArAr.age(x,exterr=exterr,i=i,sigdig=sigdig,...)
    out
}
#' @param central logical flag indicating whether each U-Th-He analysis
#'     should be considered separately (\code{central=FALSE}) or a
#'     central age should be calculated from all U-Th-He analyses
#'     together (\code{central=TRUE}).
#' @rdname age
#' @export
age.UThHe <- function(x,central=FALSE,i=NA,sigdig=NA,...){
    if (central) out <- central(x)
    else out <- UThHe.age(x,i=i,sigdig=sigdig)
    out
}
#' @rdname age
#' @export
age.fissiontracks <- function(x,central=FALSE,i=NA,sigdig=NA,exterr=TRUE,...){
    if (central) out <- central(x)
    else out <- fissiontrack.age(x,i=i,sigdig=sigdig,exterr=exterr)
    out
}
