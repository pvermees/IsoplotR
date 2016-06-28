#' Calculate isotopic ages
#'
#' Calculates U-Pb ages and propagates their analytical
#' uncertainties. Evaluates the equivalence of multiple
#' (\eqn{^{206}}Pb/\eqn{^{238}}U-\eqn{^{207}}Pb/\eqn{^{235}}U or
#' \eqn{^{207}}Pb/\eqn{^{206}}Pb-\eqn{^{206}}Pb/\eqn{^{238}}U)
#' compositions, computes the weighted mean isotopic composition and
#' the corresponding concordia age using the method of maximum
#' likelihood, computes the mswd of equivalence and concordance and
#' their respective Chi-squared p-values. Performs linear regression
#' of U-Pb data on Wetherill and Tera-Wasserburg concordia
#' diagrams. Computes the upper and lower intercept ages (for
#' Wetherill) or the lower intercept age and the
#' \eqn{^{207}}Pb/\eqn{^{206}}Pb intercept (for Tera-Wasserburg),
#' taking into account error correlations and decay constant
#' uncertainties.
#'
#' @param x a scalar containing an isotopic ratio, a two element
#'     vector containing an isotopic ratio and its standard error, or
#'     an object of class \code{UPb} or \code{detritals}.
#' @param method one of either \code{'Pb206U238'}, \code{'Pb207U235'},
#'     or \code{'Pb207Pb206'}
#' @param dcu propagate the decay constant uncertainties?
#' @param i (optional) index of a particular aliquot
#' @param ... optional arguments
#' @return if \code{x} is a scalar or a vector, returns the age using
#'     the geochronometer given by \code{method} and its standard
#'     error.
#' @rdname age
#' @export
age <- function(x,...){ UseMethod("age",x) }
#' @rdname age
#' @export
age.default <- function(x,method='Pb206U238',dcu=TRUE,...){
    if (length(x)==1) X <- c(x,0)
    else X <- x[1:2]
    if (identical(method,'Pb207U235')) out <- get.Pb207U235age(X[1],X[2],dcu)
    else if (identical(method,'Pb206U238')) out <- get.Pb206U238age(X[1],X[2],dcu)
    else if (identical(method,'Pb207Pb206')) out <- get.Pb207Pb206age(X[1],X[2],dcu)
}
#' @param concordia scalar flag indicating whether each U-Pb analysis
#'     should be considered separately (\code{concordia=1}), a
#'     concordia age should be calculated from all U-Pb analyses
#'     together (\code{concordia=2}), or a discordia line should be
#'     fit through all the U-Pb analyses (\code{concordia=2}).
#' @param wetherill boolean flag to indicate whether the data should
#'     be evaluated in Wetherill (\code{TRUE}) or Tera-Wasserburg
#'     (\code{FALSE}) space.  This option is only used when
#'     \code{concordia=2}
#' @return
#' if \code{x} has class \code{UPb} and \code{concordia=1}, returns a
#' table with the following columns: `t.75', `err[t.75]', `t.68',
#' `err[t.68]', `t.76',`err[t.76]', `t.conc', `err[t.conc]',
#' containing the 207Pb/235U-age and standard error, the
#' \eqn{^{206}}Pb/\eqn{^{238}}U-age and standard error, the
#' \eqn{^{207}}Pb/\eqn{^{206}} Pb-age and standard error, and the
#' concordia age and standard error, respectively.
#' 
#' 
#' if \code{x} has class \code{UPb} and \code{concordia=2}, returns a
#' list with the following items:
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
#' if \code{x} has class \code{UPb} and \code{concordia=3}, returns a
#' list with the following items:
#'
#' \describe{
#' \item{x}{ a two element vector with the upper and lower intercept
#' ages (if wetherill==TRUE) or the lower intercept age and
#' \eqn{^{207}}Pb/\eqn{^{206}}Pb intercept (for Tera-Wasserburg) }
#' 
#' \item{cov}{ the covariance matrix of the elements in \code{x} }
#' }
#' @examples
#' data(examples)
#' print(age(examples$UPb))
#' print(age(examples$UPb,concordia=1))
#' print(age(examples$UPb,concordia=2))
#' @rdname age
#' @export
age.UPb <- function(x,concordia=1,wetherill=TRUE,dcu=TRUE,i=NA,...){
    if (concordia==1) { out <- UPb.age(x,dcu=dcu,i=i,...) }
    else if (concordia==2) { out <- concordia.age(x,wetherill=TRUE,dcu=TRUE,...) }
    else if (concordia==3) { out <- discordia.age(x,wetherill=TRUE,dcu=TRUE) }
    out
}
#' @rdname age
#' @export
age.detritals <- function(x,...){
    x
}
#' @param isochron boolean flag indicating whether each Ar-Ar analysis
#'     should be considered separately (\code{isochron=FALSE}) or an
#'     isochron age should be calculated from all Ar-Ar analyses
#'     together (\code{isochron=TRUE}).
#' @rdname age
#' @export
age.ArAr <- function(x,isochron=FALSE,dcu=TRUE,i=NA,...){
    if (isochron){
        out <- isochron(x,plot=FALSE)
    } else {
        out <- ArAr.age(x,dcu=dcu,i=i,...)
    }
    out
}
