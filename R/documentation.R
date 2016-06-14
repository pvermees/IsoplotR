#' Example datasets for testing \code{IsoplotR}
#'
#' U-Pb and detrital zircon datasets
#'
#' \code{examples} is a list with two items
#'
#' \code{UPb}: an object of class \code{'UPb'} containing a high
#' precision U-Pb dataset packaged with Ken Ludwig's \code{Isoplot}
#' program.
#'
#' \code{DZ}: an object of class \code{'detrital'} containing a
#' detrital zircon U-Pb dataset from Namibia.
#' 
#' @name examples
#' @docType data
#' @examples
#' data(examples)
#' concordia(examples$UPb)
#' dev.new()
#' kde(examples$DZ)
#' @author Ken Ludwig and Pieter Vermeesch
#' @references
#' 
#' Ludwig, K. R. User's manual for Isoplot 3.00: a geochronological
#'     toolkit for Microsoft Excel. No. 4. Kenneth R. Ludwig, 2003.
#'
#' Vermeesch, Pieter, and Eduardo
#' Garzanti. "Making geological sense of 'Big Data' in sedimentary provenance analysis."
#' Chemical Geology 409 (2015): 20-27.
NULL
