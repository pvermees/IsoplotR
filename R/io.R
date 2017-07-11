#' Read geochronology data
#'
#' Cast a \code{.csv} file or a matrix into one of \code{IsoplotR}'s
#' data classes
#'
#' @param x either a file name (\code{.csv} format) OR a matrix
#' @param method one of \code{'U-Pb'}, \code{'Pb-Pb'}, \code{'Ar-Ar'},
#'     \code{'detritals'}, \code{Rb-Sr}, \code{Sm-Nd}, \code{Re-Os},
#'     \code{Th-U}, \code{'U-Th-He'}, \code{'fissiontracks'} or
#'     \code{'other'}
#' @param format formatting option, depends on the value of
#'     \code{method}.
#' 
#' if \code{method='U-Pb'}, then \code{format} is one of either:
#'
#' \enumerate{
#' \item{\code{7/5, s[7/5], 6/8, s[6/8], rho}}
#' \item{\code{8/6, s[8/6], 7/6, s[7/6] (, rho)}}
#' \item{\code{X=7/6, s[X], Y=7/5, s[Y], Z=6/8, s[Z] (, rho[X,Y]) (, rho[Y,Z])}}
#' }
#'
#' where optional columns are marked in round brackets
#'
#' if \code{method='Pb-Pb'}, then \code{format} is one of either:
#'
#' \enumerate{
#' \item{\code{6/4, s[6/4], 7/4, s[7/4], rho}}
#' \item{\code{4/6, s[4/6], 7/6, s[7/6], rho}}
#' \item{\code{6/4, s[6/4], 7/4, s[7/4], 7/6, s[7/6]}}
#' }
#'
#' if \code{method='Ar-Ar'}, then \code{format} is one of either:
#'
#' \enumerate{
#' \item{\code{9/6, s[9/6], 0/6, s[0/6], rho (, 39)}}
#' \item{\code{6/0, s[6/0], 9/0, s[9/0] (, rho) (, 39)}}
#' \item{\code{9/0, s[9/0], 6/0, s[6/0], 9/6, s[9/6] (, 39)}}
#' }
#'
#' if \code{method='Rb-Sr'}, then \code{format} is one of either:
#'
#' \enumerate{
#' \item{\code{Rb87/Sr86, s[Rb87/Sr86], Sr87/Sr86, s[Sr87/Sr86] (, rho)}}
#' \item{\code{Rb, s[Rb], Sr, s[Sr], Sr87/Sr86, s[Sr87/Sr86]}}
#' }
#'
#' where \code{Rb} and \code{Sr} are in ppm
#'
#' if \code{method='Sm-Nd'}, then \code{format} is one of either:
#'
#' \enumerate{
#' \item{\code{Sm147/Nd144, s[Sm147/Nd144], Nd143/Nd144, s[Nd143/Nd144] (, rho)}}
#' \item{\code{Sm, s[Sm], Nd, s[Nd], Nd143/Nd144, s[Nd143/Nd144]}}
#' }
#'
#' where \code{Sm} and \code{Nd} are in ppm
#'
#' if \code{method='Re-Os'}, then \code{format} is one of either:
#'
#' \enumerate{
#' \item{\code{Re187/Os188, s[Re187/Os188], Os187/Os188, s[Os187/Os188] (, rho)}}
#' \item{\code{Re, s[Re], Os, s[Os], Os187/Os188, s[Os187/Os188]}}
#' }
#'
#' where \code{Re} and \code{Os} are in ppm
#'
#' if \code{method='Lu-Hf'}, then \code{format} is one of either:
#'
#' \enumerate{
#' \item{\code{Lu176/Hf177, s[Lu176/Hf177], Hf176/Hf177, s[Hf176/Hf177] (, rho)}}
#' \item{\code{Lu, s[Lu], Hf, s[Hf], Hf176/Hf177, s[Hf176/Hf177]}}
#' }
#'
#' where \code{Lu} and \code{Hf} are in ppm
#'
#' if \code{method='Th-U'}, then \code{format} is one of either:
#'
#' \enumerate{
#' \item{\code{X=8/2, s[X], Y=4/2, s[Y], Z=0/2, s[Z], rho[X,Y], rho[X,Z], rho[Y,Z]}}
#' \item{\code{X=2/8, s[X], Y=4/8, s[Y], Z=0/8, s[Z], rho[X,Y], rho[X,Z], rho[Y,Z]}}
#' }
#'
#' where all values are activity ratios
#'
#' if \code{method='fissiontracks'}, then \code{format} is one of
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
#' uncertainties.}  }
#'
#' @details IsoplotR provides the following example input files:
#'
#' \itemize{
#' \item{U-Pb: \code{UPb1.csv}, \code{UPb2.csv}, \code{UPb3.csv}}
#' \item{Pb-Pb: \code{PbPb1.csv}, \code{PbPb2.csv}, \code{PbPb3.csv}}
#' \item{Ar-Ar: \code{ArAr1.csv}, \code{ArAr2.csv}, \code{ArAr3.csv}}
#' \item{Re-Os: \code{ReOs1.csv}, \code{ReOs2.csv}}
#' \item{Sm-Nd: \code{SmNd1.csv}, \code{SmNd2.csv}}
#' \item{Rb-Sr: \code{RbSr1.csv}, \code{RbSr2.csv}}
#' \item{Lu-Hf: \code{LuHf1.csv}, \code{LuHf2.csv}}
#' \item{Th-U: \code{ThU1.csv}, \code{ThU2.csv}}
#' \item{fissiontracks: \code{FT1.csv}, \code{FT2.csv}, \code{FT3.csv}}
#' \item{U-Th-He: \code{UThHe.csv}, \code{UThSmHe.csv}}
#' \item{detritals: \code{Namib.csv}}
#' \item{other: \code{MountTom.csv}, \code{average.csv}, \code{spectrum.csv}}
#' }
#' 
#' The contents of these files can be viewed using the
#' \code{system.file(...)} function.
#'
#' @param ... optional arguments to the \code{read.csv} function
#' @return an object of class \code{UPb}, \code{PbPb}, \code{ArAr},
#'     \code{UThHe}, \code{ReOs}, \code{SmNd}, \code{RbSr},
#'     \code{LuHf}, \code{detritals}, \code{fissiontracks} or
#'     \code{other}
#' @examples
#' file.show(system.file("spectrum.csv",package="IsoplotR"))
#'
#' f1 <- system.file("UPb1.csv",package="IsoplotR")
#' d1 <- read.data(f1,method="U-Pb",format=1)
#' concordia(d1)
#'
#' f2 <- system.file("ArAr1.csv",package="IsoplotR")
#' d2 <- read.data(f2,method="Ar-Ar",format=1)
#' agespectrum(d2)
#'
#' f3 <- system.file("ReOs1.csv",package="IsoplotR")
#' d3 <- read.data(f3,method="Re-Os",format=1)
#' isochron(d2)
#'
#' f4 <- system.file("FT1.csv",package="IsoplotR")
#' d4 <- read.data(f4,method="fissiontracks",format=1)
#' radialplot(d4)
#'
#' f5 <- system.file("UThSmHe.csv",package="IsoplotR")
#' d5 <- read.data(f5,method="U-Th-He")
#' helioplot(d5)
#'
#' f6 <- system.file("ThU2.csv",package="IsoplotR")
#' d6 <- read.data(f6,method="Th-U",format=2)
#' evolution(d6)
#'
#' #  one detrital zircon U-Pb file (detritals.csv)
#' f7 <- system.file("Namib.csv",package="IsoplotR")
#' d7 <- read.data(f7,method="detritals")
#' kde(d7)
#'
#' #  three 'other' files (MountTom.csv, spectrum.csv, average.csv)
#' f8 <- system.file("MountTom.csv",package="IsoplotR")
#' d8 <- read.data(f8,method="other")
#' radialplot(d8)
#'
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
    } else if (identical(method,'Pb-Pb')){
        out <- as.PbPb(x,format)
    } else if (identical(method,'Ar-Ar')){
        out <- as.ArAr(x,format)
    } else if (identical(method,'Re-Os')){
        out <- as.ReOs(x,format)
    } else if (identical(method,'Rb-Sr')){
        out <- as.RbSr(x,format)
    } else if (identical(method,'Sm-Nd')){
        out <- as.SmNd(x,format)
    } else if (identical(method,'Lu-Hf')){
        out <- as.LuHf(x,format)
    } else if (identical(method,'Th-U')){
        out <- as.ThU(x,format)
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

as.UPb <- function(x,format=3){
    out <- list()
    class(out) <- "UPb"
    out$x <- NA
    out$format <- format
    nc <- ncol(x)
    nr <- nrow(x)
    X <- shiny2matrix(x,2,nr,nc)
    cnames <- NULL
    if (format == 1 & nc == 5){
        cnames <- c('Pb207U235','errPb207U235',
                    'Pb206U238','errPb206U238','rhoXY')
        out$x <- X
    } else if (format == 2 & nc %in% c(4,5)) {
        cnames <- c('U238Pb206','errU238Pb206',
                    'Pb207Pb206','errPb207Pb206','rhoXY')
        if (nc == 4){
            rho <- rep(0,nr-1)
            out$x <- cbind(X,rho)
        } else {
            i <- which(is.na(X[,5]))
            X[i,5] <- 0
            out$x <- X
        }
    } else if (format == 3 & nc %in% c(6,7,8)){
        cnames <- c('Pb207U235','errPb207U235',
                    'Pb206U238','errPb206U238',
                    'Pb207Pb206','errPb207Pb206',
                    'rhoXY','rhoYZ')
        if (nc == 8){
            rhoXY <- X[,7]
            rhoYZ <- X[,8]
            i <- which(is.na(rhoXY))
            j <- which(is.na(rhoYZ))
        } else if (nc == 7){
            rhoXY <- X[,7]
            rhoYZ <- rep(0,nr-1)
            i <- which(is.na(rhoXY))
            j <- 1:(nr-1)
            X <- cbind(X,rhoYZ)
        } else {
            rhoXY <- rep(0,nr-1)
            rhoYZ <- rep(0,nr-1)
            i <- 1:(nr-1)
            j <- 1:(nr-1)
            X <- cbind(X,rhoXY,rhoYZ)
        }
        X[i,7] <- get.cor.75.68(X[i,1],X[i,2],X[i,3],X[i,4],X[,5],X[i,6])
        X[j,8] <- get.cor.68.76(X[j,1],X[j,2],X[j,3],X[j,4],X[,5],X[j,6])
        out$x <- X
    }
    colnames(out$x) <- cnames
    out
}
get.cor.75.68 <- function(Pb207U235,errPb207U235,Pb206U238,errPb206U238,
                          Pb207Pb206,errPb207Pb206){
    get.cor.div(Pb207U235,errPb207U235,
                Pb206U238,errPb206U238,
                Pb207Pb206,errPb207Pb206)
}
get.cor.68.76 <- function(Pb207U235,errPb207U235,Pb206U238,errPb206U238,
                          Pb207Pb206,errPb207Pb206){
    get.cor.mult(Pb206U238,errPb206U238,
                 Pb207Pb206,errPb207Pb206,
                 Pb207U235,errPb207U235)
}
as.PbPb <- function(x,format=1){
    out <- list()
    class(out) <- "PbPb"
    out$x <- NA
    out$format <- format
    nc <- ncol(x)
    nr <- nrow(x)
    X <- shiny2matrix(x,2,nr,nc)
    cnames <- NULL
    if (format == 1 & nc == 5){
        cnames <- c('Pb206Pb204','errPb206Pb204',
                    'Pb207Pb204','errPb207Pb204','rho')
    } else if (format == 2 & nc == 5) {
        cnames <- c('Pb204Pb206','errPb204Pb206',
                    'Pb207Pb206','errPb207Pb206','rho')
    } else if (format == 3 & nc == 6){
        cnames <- c('Pb206Pb204','errPb206Pb204',
                    'Pb207Pb204','errPb207Pb204',
                    'Pb207Pb206','errPb207Pb206')
    }
    out$x <- X
    colnames(out$x) <- cnames
    out
}
as.ArAr <- function(x,format=3){
    out <- list()
    class(out) <- "ArAr"
    out$format <- format
    out$J <- as.numeric(x[2,1:2])
    nc <- ncol(x)
    nr <- nrow(x)
    bi <- 4 # begin index
    X <- shiny2matrix(x,bi,nr,nc)
    if (format==3 & nc %in% c(6,7)){
        if (nc==7){
            out$x <- X
        } else {
            ns <- nr-bi+1 # number of samples
            Ar39 <- rep(1/ns,ns)
            out$x <- cbind(X,Ar39)            
        }
        colnames(out$x) <- c('Ar39Ar40','errAr39Ar40',
                             'Ar36Ar40','errAr36Ar40',
                             'Ar39Ar36','errAr39Ar36','Ar39')
    } else if (nc %in% c(4,5,6)){
        if (nc==6){
            out$x <- X
        } else {
            ns <- nr-bi+1 # number of samples
            Ar39 <- rep(1/ns,ns)
            out$x <- cbind(X,Ar39)
        }
        if (nc==4) {
            rho <- 0
        }
        if (format==1) {
            colnames(out$x) <- c('Ar39Ar36','errAr39Ar36',
                                 'Ar40Ar36','errAr40Ar36',
                                 'rho','Ar39')
        } else {
            colnames(out$x) <- c('Ar39Ar40','errAr39Ar40',
                                 'Ar36Ar40','errAr36Ar40',
                                 'rho','Ar39')
        }
    }
    out
}
as.RbSr <- function(x,format=2){
    cnames1 <- c('Rb87Sr86','errRb87Sr86',
                 'Sr87Sr86','errSr87Sr86','rho')
    cnames2 <- c('Rbppm','errRbppm','Srppm','errSrppm',
                 'Sr87Sr86','errSr87Sr86')
    as.PD(x,"RbSr",cnames1,cnames2,format)
}
as.ReOs <- function(x,format=2){
    cnames1 <- c('Re187Os188','errRe187Os188',
                 'Os187Os188','errOs187Os188','rho')
    cnames2 <- c('Reppm','errReppm','Osppm','errOsppm',
                 'Os187Os188','errOs187Os188')
    as.PD(x,"ReOs",cnames1,cnames2,format)
}
as.SmNd <- function(x,format=2){
    cnames1 <- c('Sm143Nd144','errSm143Nd144',
                 'Nd143Nd144','errNd143Nd144','rho')
    cnames2 <- c('Smppm','errSmppm','Ndppm','errNdppm',
                 'Nd143Nd144','errNd143Nd144')
    as.PD(x,"SmNd",cnames1,cnames2,format)
}
as.LuHf <- function(x,format=2){
    cnames1 <- c('Lu176Hf177','errLu176Hf177',
                 'Hf176Hf177','errHf176Hf177','rho')
    cnames2 <- c('Luppm','errLuppm','Hfppm','errHfppm',
                 'Hf176Hf177','errHf176Hf177')
    as.PD(x,"LuHf",cnames1,cnames2,format)
}
as.PD <- function(x,classname,colnames1,colnames2,format){
    out <- list()
    class(out) <- c(classname,'PD')
    out$x <- NA
    out$format <- format
    nc <- ncol(x)
    nr <- nrow(x)
    X <- shiny2matrix(x,2,nr,nc)
    if (format == 1 & nc %in% c(4,5)){
        colnames(X) <- colnames1
        if (nc == 4){
            rho <- rep(0,nr-1)
            out$x <- cbind(X,rho)
        } else if (nc == 5){
            i <- which(is.na(X[,5]))
            X[i,5] <- 0
            out$x <- X
        }
    } else if (format == 2 & nc == 6){
        colnames(X) <- colnames2
        out$x <- X
    }
    out
}
as.ThU <- function(x,format=1){
    out <- list()
    class(out) <- "ThU"
    out$x <- NA
    out$format <- format
    nc <- ncol(x)
    nr <- nrow(x)
    X <- shiny2matrix(x,2,nr,nc)
    cnames <- NULL
    if (format == 1 & nc == 9){
        cnames <- c('U238Th232','errU238Th232',
                    'U234Th232','errU234Th232',
                    'Th230Th232','errTh230Th232',
                    'rhoXY','rhoXZ','rhoYZ')
    } else if (format == 2 & nc == 9) {
        cnames <- c('Th232U238','errTh232U238',
                    'U234U238','errU234U238',
                    'Th230U238','errTh230U238',
                    'rhoXY','rhoXZ','rhoYZ')
    }
    out$x <- X
    colnames(out$x) <- cnames
    out
}
as.UThHe <- function(x){
    nc <- ncol(x)
    nr <- nrow(x)
    out <- matrix(0,nr-1,nc)
    out[1:(nr-1),1:nc] <- shiny2matrix(x,2,nr,nc)
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
        X <- shiny2matrix(x,6,nr,2)
        out$x <- shiny2matrix(x,6,nr,2)
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
    X <- shiny2matrix(x,2,nr,nc)
    colnames(X) <- snames
    for (sname in snames){
        out[[sname]] = X[!is.na(X[,sname]),sname]
    }
    out
}
as.other <- function(x){
    nc <- ncol(x)
    has.header <- is.na(suppressWarnings(as.numeric(x[1,1])))
    if (has.header) x <- x[-1,]
    matrix(as.numeric(x),ncol=nc)
}

# x = a numerical vector, br = length of the preamble with parameters
# nr = number of rows, nc = number of columns
shiny2matrix <- function(x,br,nr,nc){
    suppressWarnings(
        return(matrix(as.numeric(x[(br:nr),]),nr-br+1,nc))
    )
}
