#' Read geochronology data
#'
#' Cast a \code{.csv} file or a matrix into one of \code{IsoplotR}'s
#' data classes
#'
#' @param x either a file name (\code{.csv} format) OR a matrix
#' @param method one of \code{'U-Pb'}, \code{'Ar-Ar'},
#'     \code{'detritals'} \code{'U-Th-He'}, \code{'fissiontracks'} or
#'     \code{'other'}
#' @param format formatting option, depends on the value of
#'     \code{method}.
#' 
#' if \code{method = 'Ar-Ar'}, then \code{format} is one of either:
#'
#' \enumerate{
#' \item{\code{39/40, s[39/40], 36/40, s[36/40], 39/36, s[39/36]}}
#' \item{\code{39, 39/40, s[39/40], 36/40, s[36/40], 39/36, s[39/36]}}
#' }
#'
#' if \code{method = 'fissiontracks'}, then \code{format} is one of
#' either:
#'
#' \enumerate{
#' \item{the External Detector Method (EDM), which requires a
#' \eqn{\zeta}-calibration constant and its uncertainty, the induced
#' track density in a dosimeter glass, and a table with the
#' spontaneous and induced track densities.}
#'
#' \item{LA-ICP-MS-based fission track data using the
#' \eqn{\zeta}-calibration method, which requires a 'session
#' \eqn{\zeta}' and its uncertainty and a table with the number of
#' spontaneous tracks, the area over which these were counted and one
#' or more U/Ca- or U-concentration measurements and their analytical
#' uncertainties.}
#'
#' \item{LA-ICP-MS-based fission track data using the 'absolute
#' dating' method, which only requires a table with the the number of
#' spontaneous tracks, the area over which these were counted and one
#' or more U/Ca- or U-concentration measurements and their analytical
#' uncertainties.}
#' }
#'
#' @details
#' Example input files can be found by using \code{R}'s
#' \code{system.file(...)} function:
#'
#' \enumerate{
#' \item \code{method = 'U-Pb'}:
#'
#' \code{file.show(system.file("UPb.csv",package="IsoplotR"))}
#'
#' \item \code{method = 'Ar-Ar'} and \code{format = 1}:
#'
#' \code{file.show(system.file("ArAr1.csv",package="IsoplotR"))}
#'
#' \item \code{method = 'Ar-Ar'} and \code{format = 2}:
#'
#' \code{file.show(system.file("ArAr2.csv",package="IsoplotR"))}
#'
#' \item \code{method = 'Re-Os'}:
#'
#' \code{file.show(system.file("ReOs.csv",package="IsoplotR"))}
#'
#' \item \code{method = 'U-Th-He'}:
#'
#' \code{file.show(system.file("UThHe.csv",package="IsoplotR"))}
#'
#' \item \code{method = 'fissiontracks'} and \code{format = 1}:
#'
#' \code{file.show(system.file("FT1.csv",package="IsoplotR"))}
#' 
#' \item \code{method = 'fissiontracks'} and \code{format = 2}:
#'
#' \code{file.show(system.file("FT2.csv",package="IsoplotR"))}
#' 
#' \item \code{method = 'fissiontracks'} and \code{format = 3}:
#'
#' \code{file.show(system.file("FT3.csv",package="IsoplotR"))}
#' 
#' \item \code{method = 'detritals'}:
#'
#' \code{file.show(system.file("DZ.csv",package="IsoplotR"))}
#' }
#' @param ... optional arguments to the \code{read.csv} function
#' @return an object of class \code{UPb}, \code{ArAr}, \code{UThHe},
#'     \code{detritals}, \code{fissiontracks} or \code{other}
#' @examples
#' # load one of the built-in .csv files:
#' data(examples)
#' concordia(examples$UPb)
#' @rdname read.data
#' @export
read.data <- function(x,...){ UseMethod("read.data",x) }
#' @rdname read.data
#' @export
read.data.default <- function(x,method='U-Pb',format=1,...){
    X <- as.matrix(utils::read.table(x,sep=',',...))
    read.data.matrix(X,method=method,format=format)
}
#' @rdname read.data
#' @export
read.data.matrix <- function(x,method='U-Pb',format=1,...){
    if (identical(method,'U-Pb')){
        out <- as.UPb(x,format)
    } else if (identical(method,'Ar-Ar')){
        out <- as.ArAr(x,format)
    } else if (identical(method,'Re-Os')){
        out <- as.ReOs(x,format)
    } else if (identical(method,'U-Th-He')){
        out <- as.UThHe(x)
    } else if (identical(method,'fissiontracks')){
        out <- as.fissiontracks(x,format)
    } else if (identical(method,'detritals')){
        out <- as.detritals(x)
    } else if (identical(method,'other')){
        out <- as.other(x)
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
        X <- x2X(x,2,nr,nc)
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
as.ArAr <- function(x,format=2){
    out <- list()
    class(out) <- "ArAr"
    out$format <- format
    out$J <- as.numeric(x[2,1:2])
    nc <- ncol(x)
    nr <- nrow(x)
    if (format == 1 & nc == 6){
        X <- x2X(x,4,nr,nc)
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
    } else if (format == 2 & nc == 7){
        X <- x2X(x,4,nr,nc)
        colnames(X) <- c('Ar39',
                         'Ar39Ar40','errAr39Ar40',
                         'Ar36Ar40','errAr36Ar40',
                         'Ar39Ar36','errAr39Ar36')
        Ar40Ar36 <- 1/X[,'Ar36Ar40']
        errAr40Ar36 <- X[,'errAr36Ar40']/X[,'Ar36Ar40']^2
        out$x <- cbind(X,Ar40Ar36,errAr40Ar36)
        colnames(out$x) <- c('Ar39',
                             'Ar39Ar40','errAr39Ar40',
                             'Ar36Ar40','errAr36Ar40',
                             'Ar39Ar36','errAr39Ar36',
                             'Ar40Ar36','errAr40Ar36')
    }
    out
}
as.ReOs <- function(x,format=1){
    out <- list()
    class(out) <- "ReOs"
    out$x <- NA
    out$format <- format
    nc <- ncol(x)
    nr <- nrow(x)
    if (format == 1 & nc == 6){
        X <- x2X(x,2,nr,nc)
        colnames(X) <- c('Reppm','errReppm',
                         'Osppt','errOsppt',
                         'Os187Os188','errOs187Os188')
        out$x <- X
    }
    out
}
as.UThHe <- function(x){
    nc <- ncol(x)
    nr <- nrow(x)
    out <- matrix(0,nr-1,nc)
    out[1:(nr-1),1:nc] <- x2X(x,2,nr,nc)
    if (nc==8) {
        colnames(out) <- c('He','errHe','U','errU','Th','errTh','Sm','errSm')
    } else if (nc==6) {
        colnames(out) <- c('He','errHe','U','errU','Th','errTh')
    } else {
        stop("Input table must have 6 or 8 columns.")
    }
    class(out) <- append(class(out),"UThHe")
    out
}
as.fissiontracks <- function(x,format=1){
    nr <- nrow(x)
    nc <- ncol(x)
    out <- list()
    class(out) <- "fissiontracks"
    out$format <- format
    if (format==1){
        out$zeta <- as.numeric(x[2,1:2])
        out$rhoD <- as.numeric(x[4,1:2])
        X <- x2X(x,6,nr,2)
        out$x <- x2X(x,6,nr,2)
        colnames(out$x) <- c('Ns','Ni')
    } else {
        if (format==2){
            out$zeta <- as.numeric(x[2,1:2])
            out$spotSize <- as.numeric(x[4,1])
            si <- 6 # start index
        } else {
            out$spotSize <- as.numeric(x[2,1])
            si <- 4
        }
        ns <- nr-si+1
        Ns <- as.numeric(x[si:nr,1])
        A <- as.numeric(x[si:nr,2])
        out$Ns <- Ns
        out$A <- A
        out$U <- list()
        out$sU <- list()
        j <- seq(from=3,to=nc-1,by=2)
        for (i in 1:ns){
            U <- as.numeric(x[i+si-1,j])
            out$U[[i]] <- U
            sU <- as.numeric(x[i+si-1,j+1])
            out$sU[[i]] <- sU
        }
    }
    out
}
as.detritals <- function(x){
    out <- list()
    class(out) <- "detritals"
    nr <- nrow(x)
    nc <- ncol(x)
    snames <- x[1,]
    X <- x2X(x,2,nr,nc)
    colnames(X) <- snames
    for (sname in snames){
        out[[sname]] = X[!is.na(X[,sname]),sname]
    }
    out
}
as.RbSr <- function(x){

}
as.other <- function(x){
    nc <- ncol(x)
    has.header <- is.na(suppressWarnings(as.numeric(x[1,1])))
    if (has.header) x <- x[-1,]
    matrix(as.numeric(x),ncol=nc)
}

x2X <- function(x,br,nr,nc){
    suppressWarnings(
        return(matrix(as.numeric(x[(br:nr),]),nr-br+1,nc))
    )
}
