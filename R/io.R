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
#' 
#' \code{1}: 7/6, s[7/6], 6/8, s[6/8], 7/5, s[7/5]
#' 
#' (other formats will be added later)
#' @param ... optional arguments to the \code{read.csv} function
#' @return an object of class \code{'UPb'}, \code{'ArAr'},
#'     \code{'RbSr'}, \code{'SmNd'}, \code{'ReOs'}, \code{'UThHe'},
#'     \code{'fission'}, \code{'cosmogenics'}, or \code{'other'}
#' @examples
#' # load one of the built-in .csv files:
#' data(examples)#fname <- system.file("UPb.csv",package="IsoplotR")
#' #UPb <- read.data(fname,'U-Pb')
#' concordia(examples$UPb)
#' @rdname read.data
#' @export
read.data <- function(x,...){ UseMethod("read.data",x) }
#' @rdname read.data
#' @export
read.data.default <- function(x,method='Pb206U238',format=1,...){
    X <- as.matrix(utils::read.table(x,sep=',',...))
    read.data.matrix(X,method=method,format=format)
}
#' @rdname read.data
#' @export
read.data.matrix <- function(x,method='U-Pb',format=1,...){
    if (identical(method,'U-Pb')){
        out <- as.UPb(x,format)
    } else if (identical(method,'Ar-Ar')){
        out <- as.ArAr(x)
    } else if (identical(method,'detritals')){
        out <- as.detritals(x)
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
        X <- matrix(as.numeric(x[(2:nr),]),nr-1,nc)
        colnames(X) <- c('Pb207Pb206','errPb207Pb206',
                         'Pb206U238','errPb206U238',
                         'Pb207U235','errPb207U235')
        U238Pb206 <- 1/X[,'Pb206U238']
        errU238Pb206 <- X[,'errPb206U238']/X[,'Pb206U238']^2
        out$x <- cbind(X[,c('Pb207U235','errPb207U235','Pb206U238','errPb206U238')],
                       U238Pb206, errU238Pb206,
                       X[,c('Pb207Pb206','errPb207Pb206')])
        colnames(out$x) <- c('Pb207U235','errPb207U235',
                             'Pb206U238','errPb206U238',
                             'U238Pb206','errU238Pb206',
                             'Pb207Pb206','errPb207Pb206')
    }
    out
}
as.ArAr <- function(x,format=1){
    out <- list()
    class(out) <- "ArAr"
    out$format <- format
    out$J <- as.numeric(x[2,1:2])
    nc <- ncol(x)
    nr <- nrow(x)
    if (format == 1 & nc == 6){
        X <- matrix(as.numeric(x[4:nr,]),nr-3,nc)
        colnames(X) <- c('Ar39Ar40','errAr39Ar40',
                         'Ar36Ar40','errAr36Ar40',
                         'Ar39Ar36','errAr39Ar36')
        Ar40Ar36 <- 1/X[,'Ar36Ar40']
        errAr40Ar36 <- X[,'errAr36Ar40']/X[,'Ar36Ar40']^2
        out$x <- cbind(X,Ar40Ar36,errAr40Ar36)
        colnames(out$x) <- c('Ar39Ar40','errAr39Ar40',
                             'Ar36Ar40','errAr36Ar40',
                             'Ar39Ar36','errAr39Ar36',
                             'Ar40Ar36','errAr40Ar36')
    }
    out
}
as.detritals <- function(x){
    nr <- nrow(x)
    nc <- ncol(x)
    snames <- x[1,]
    X <- matrix(as.numeric(x[(2:nr),]),nr-1,nc)
    colnames(X) <- snames
    out <- list()
    class(out) <- "detritals"
    for (sname in snames){
        out[[sname]] = X[!is.na(X[,sname]),sname]
    }
    out
}
