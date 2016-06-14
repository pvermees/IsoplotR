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
age.default <- function(x,...){ return(x) }
#' @param concordia one of either
#'
#' 0: consider each U-Pb analysis separately
#'
#' 1: calculate a concordia age from all U-Pb analyses together
#'
#' 2: fit a discordia line through all the U-Pb analyses
#' @return if \code{x} has class \code{UPb}, \code{concordia}=0,
#'     returns a table with the following columns: 't.75',
#'     'err[t.75]','t.68','err[t.68]', 't.76','err[t.76]',
#'     't.conc','err[t.conc]', containing the 207Pb/235U-age
#'     and standard error, the 206Pb/238U-age and standard error,
#'     the 207Pb/206Pb-age and standard error, and the concordia
#'     age and standard error, respectively.
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
