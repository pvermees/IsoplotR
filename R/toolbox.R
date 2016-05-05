#' Get the covariance matrix of a sample
#'
#' Returns the covariance matrix of the i\eqn{^{th}} sample
#'
#' @param x an object of class \code{UPb}
#' @param i the index of the sample of interest
#' @return a covariance matrix of size [3x3]
#' @examples
#' data(UPb)
#' get.covmat(UPb,2)
#' @export
get.covmat <- function(x,i){
    if (methods::is(x,'UPb')){
        return(get.covmat.UPb(x,i))
    }
}

get.ages <- function(x){
    if (methods::is(x,'UPb')){
        get.ages.UPb(x)
    }
}

#' Calculate the isotopic ratio for a given age
#'
#' Predict the daughter/parent ratio for a given U-Pb age
#'
#' @param age the geological age [Ma]
#' @param method currently only \code{'U-Pb'}, other chronometers will
#'     be added later
#' @return a two element list containing:
#'
#' \code{x}: a vector with predicted isotopic ratios
#'
#' \code{cov}: the covariance matrix of the predicted isotopic ratios
#' taking into account decay constant uncertainties
#' @examples
#' get.ratios(4567,'U-Pb')
#' @export
get.ratios <- function(age,method='U-Pb'){
    if (identical(method,'U-Pb')){
        return(get.ratios.UPb(age))
    }
}
