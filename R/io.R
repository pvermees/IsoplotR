#' Read geochronology data
#'
#' Cast a .csv file into one of \code{IsoplotR}'s data classes
#'
#' @param fname file name (.csv format)
#' @param method one of \code{'U-Pb'}, \code{'Ar-Ar'}, \code{'Rb-Sr'},
#'     \code{'Sm-Nd'}, \code{'Re-Os'}, \code{'U-Th-He'},
#'     \code{'fission tracks'}, \code{'cosmogenic nuclides'} or
#'     \code{'other'}
#' @param format formatting option, depends on the value of
#'     \code{method}. If \code{method = 'U-Pb'}, then \code{format} is
#'     one of either:
#'
#' \code{1}: 7/6, s[7/6], 6/8, s[6/8], 7/5, s[7/5]
#' @param ... optional arguments to the \code{read.csv} function
#' @return an object of class \code{'UPb'}, \code{'ArAr'},
#'     \code{'RbSr'}, \code{'SmNd'}, \code{'ReOs'}, \code{'UThHe'},
#'     \code{'fission'}, \code{'cosmogenics'}, or \code{'other'}
#' @examples
#' # load one of the built-in .csv files:
#' fname <- system.file("UPb.csv",package="IsoplotR")
#' UPb <- read.data(fname,'U-Pb')
#' concordia.plot(UPb)
#' @export
read.data <- function(fname,method='U-Pb',format=1,...){
    x <- utils::read.csv(fname,...)
    read.matrix(x,method,format)
}

#' Read geochronology data
#'
#' Cast a matrix into one of \code{IsoplotR}'s data classes
#'
#' @param x a matrix
#' @param method see \code{read.data} for details
#' @param format see \code{read.data} for details
#' @return see \code{read.data} for details
#' @examples
#' # load one of the built-in .csv files:
#' fname <- system.file("UPb.csv",package="IsoplotR")
#' dat <- read.csv(fname,header=TRUE)
#' UPb <- read.matrix(dat,method='U-Pb',format=1)
#' concordia.plot(UPb)
#' @export
read.matrix <- function(x,method='U-Pb',format=1){
    if (identical(method,'U-Pb')){
        return(as.UPb(x,format))
    }
}

as.UPb <- function(x,format=1){
    out <- list()
    out$x <- NA
    out$format <- format
    nc <- ncol(x)
    nr <- nrow(x)
    if (format == 1 & nc == 6){
        out$x <- as.matrix(x)
        colnames(out$x) <- c('Pb207Pb206','errPb207Pb206',
                             'Pb206U238','errPb206U238',
                             'Pb207U235','errPb207U235')
    }
    class(out) <- "UPb"
    out
}
