#' Calculate isotopic ages
#'
#' Calculates U-Pb and Ar-Ar ages
#'
#' @param x an object of class \code{UPb}
#' @param i (optional) index of a particular aliquot
#' @param ... optional arguments
#' @rdname age
#' @export
age <- function(x,...){ UseMethod("age",x) }
#' @rdname age
#' @export
age.default <- function(x,...){ stop('invalid input') }
#' @param concordia one of either
#'
#' 0: consider each U-Pb analysis separately
#'
#' 1: calculate a concordia age from all U-Pb analyses together
#'
#' 2: fit a discordia line through all the U-Pb analyses
#' @return if \code{x} has class \code{UPb}, \code{concordia}=0, and
#'     \code{i=NA}, returns a table with the following columns:
#'     '7/5-age', 'err[7/5-age]', '6/8-age', 'err[6/8-age]',
#'     '7/6-age', 'err[7/6-age]', 'concordia-age', and
#'     'err[concordia-age]'.
#'
#' if \code{x} has class \code{UPb}, \code{concordia}=0, and
#'     \code{i!=NA}, returns a list with two items (\code{x} and
#'     \code{cov}) containing the '7/5-age', '6/8-age', '7/6-age', and
#'     'concordia-age', and their covariance matrix, respectively.
#'
#' if \code{x} has class \code{UPb} and \code{concordia=1}, returns
#' the output of \link{concordia.age}
#'
#' if \code{x} has class \code{UPb} and \code{concordia=2}, returns
#' the output of \link{discordia.age}
#' @examples
#' data(examples)
#' age(examples$UPb)
#' @rdname age
#' @export
age.UPb <- function(x,concordia=0,i=NA,...){
    if (concordia==0) { UPb.age(x,i=i,...) }
    else if (concordia==1) { concordia.age(x,...) }
    else if (concordia==2) { discordia.age(x,...) }
}
