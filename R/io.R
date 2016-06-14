#' Read geochronology data
#'
#' Cast a \code{.csv} file or a matrix into one of \code{IsoplotR}'s
#' data classes
#'
#' @param x a file name (\code{.csv} format) or matrix
#' @param method one of \code{'U-Pb'}, \code{'Ar-Ar'}, \code{'Rb-Sr'},
#'     \code{'Sm-Nd'}, \code{'Re-Os'}, \code{'U-Th-He'},
#'     \code{'fission tracks'}, \code{'cosmogenic nuclides'} or
#'     \code{'other'}
#' @param format formatting option, depends on the value of
#'     \code{method}. If \code{method = 'U-Pb'}, then \code{format} is
#'     one of either:
#' \code{1}: 7/6, s[7/6], 6/8, s[6/8], 7/5, s[7/5]
#' (other formats will be added later)
#' @param ... optional arguments to the \code{read.csv} function
#' @return an object of class \code{'UPb'}, \code{'ArAr'},
#'     \code{'RbSr'}, \code{'SmNd'}, \code{'ReOs'}, \code{'UThHe'},
#'     \code{'fission'}, \code{'cosmogenics'}, or \code{'other'}
#' @examples
#' # load one of the built-in .csv files:
#' fname <- system.file("UPb.csv",package="IsoplotR")
#' UPb <- read.data(fname,'U-Pb')
#' concordia(UPb)
#' @rdname read.data
#' @export
read.data <- function(x,...){ UseMethod("read.data",x) }
#' @rdname read.data
#' @export
read.data.default <- function(x,method='Pb206U238',format=1,...){
    X <- utils::read.csv(x,...)
    read.data.matrix(X,method=method,format=format)
}
#' @rdname read.data
#' @export
read.data.matrix <- function(x,method='U-Pb',format=1,...){
    if (identical(method,'U-Pb')){
        out <- as.UPb(x,format)
    } else if (identical(method,'detritals')){
        out <- as.detritals(x)
    }
    out    
}
#' @rdname read.data
#' @export
read.data.data.frame <- function(x,method='U-Pb',format=1,...){
    read.data.matrix(as.matrix(x),method=method,format=format,...)
}

as.detritals <- function(x){
    snames <- colnames(x)
    out <- list()
    class(out) <- "detritals"
    for (sname in snames){
        out[[sname]] = x[!is.na(x[,sname]),sname]
    }
    out
}

as.UPb <- function(x,format=1){
    out <- list()
    class(out) <- "UPb"
    out$x <- NA
    out$format <- format
    nc <- ncol(x)
    nr <- nrow(x)
    if (format == 1 & nc == 6){
        X <- as.matrix(x)
        colnames(X) <- c('Pb207Pb206','errPb207Pb206',
                         'Pb206U238','errPb206U238',
                         'Pb207U235','errPb207U235')
        out$x <- cbind(X[,c('Pb207U235','errPb207U235','Pb206U238','errPb206U238')],
                       1/X[,'Pb206U238'], sqrt(X[,'errPb206U238'])/X[,'Pb206U238'],
                       X[,c('Pb207Pb206','errPb207Pb206')])
        colnames(out$x) <- c('Pb207U235','errPb207U235',
                             'Pb206U238','errPb206U238',
                             'U238Pb206','errU238Pb206',
                             'Pb207Pb206','errPb207Pb206')
    }
    out
}
