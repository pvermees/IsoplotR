#' Read geochronology data
#'
#' Cast a \code{.csv} file or a matrix into one of \code{IsoplotR}'s
#' data classes
#'
#' @details IsoplotR provides the following example input files:
#'
#' \itemize{
#' \item{U-Pb: \code{UPb1.csv}, \code{UPb2.csv}, \code{UPb3.csv},
#' \code{UPb4.csv}, \code{UPb5.csv}, \code{UPb6.csv} }
#' \item{Pb-Pb: \code{PbPb1.csv}, \code{PbPb2.csv}, \code{PbPb3.csv} }
#' \item{Ar-Ar: \code{ArAr1.csv}, \code{ArAr2.csv}, \code{ArAr3.csv}}
#' \item{K-Ca: \code{KCa1.csv}, \code{KCa2.csv}}, 
#' \item{Re-Os: \code{ReOs1.csv}, \code{ReOs2.csv}}
#' \item{Sm-Nd: \code{SmNd1.csv}, \code{SmNd2.csv}}
#' \item{Rb-Sr: \code{RbSr1.csv}, \code{RbSr2.csv}}
#' \item{Lu-Hf: \code{LuHf1.csv}, \code{LuHf2.csv}}
#' \item{Th-U: \code{ThU1.csv}, \code{ThU2.csv}, \code{ThU3.csv},
#' \code{ThU4.csv}}
#' \item{fissiontracks: \code{FT1.csv}, \code{FT2.csv},
#' \code{FT3.csv}}
#' \item{U-Th-He: \code{UThHe.csv}, \code{UThSmHe.csv}}
#' \item{detritals: \code{DZ.csv}}
#' \item{other: \code{LudwigMixture.csv}, \code{LudwigMean.csv},
#' \code{LudwigKDE.csv} \code{LudwigSpectrum.csv}}
#' }
#'
#' The contents of these files can be viewed using the
#' \code{system.file(...)} function. For example, to read the
#' \code{ArAr1.csv} file:
#'
#' \code{fname <- system.file('ArAr1.csv',package='IsoplotR')}
#'
#' \code{ArAr <- read.data(fname,method='Ar-Ar',format=1)}
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
#' \item{\code{7/5, err[7/5], 6/8, err[6/8], rho}}
#' \item{\code{8/6, err[8/6], 7/6, err[7/6] (, rho)}}
#' \item{\code{X=7/6, err[X], Y=6/8, err[Y], Z=7/6, err[Z]}
#' \code{(, rho[X,Y]) (, rho[Y,Z])}}
#' \item{\code{X=7/5, err[X], Y=6/8, err[Y], Z=4/8, }
#' \code{rho[X,Y], rho[X,Z], rho[Y,Z]}}
#' \item{\code{X=8/6, err[X], Y=7/6, err[Y], Z=4/6, }
#' \code{rho[X,Y], rho[X,Z], rho[Y,Z]}}
#' \item{\code{7/5, err[7/5], 6/8, err[6/8], 4/8, err[4/8], }
#' \code{7/6, err[7/6], 4/7, err[4/7], 4/6, err[4/6]}}
#' }
#'
#' where optional columns are marked in round brackets
#'
#' if \code{method='Pb-Pb'}, then \code{format} is one of either:
#'
#' \enumerate{
#' \item{\code{6/4, err[6/4], 7/4, err[7/4], rho}}
#' \item{\code{4/6, err[4/6], 7/6, err[7/6], rho}}
#' \item{\code{6/4, err[6/4], 7/4, err[7/4], 7/6, err[7/6]}}
#' }
#'
#' if \code{method='Ar-Ar'}, then \code{format} is one of either:
#'
#' \enumerate{
#' \item{\code{9/6, err[9/6], 0/6, err[0/6], rho (, 39)}}
#' \item{\code{6/0, err[6/0], 9/0, err[9/0] (, rho) (, 39)}}
#' \item{\code{9/0, err[9/0], 6/0, err[6/0], 9/6, err[9/6] (, 39)}}
#' }
#'
#' if \code{method='K-Ca'}, then \code{format} is one of either:
#'
#' \enumerate{
#' \item{\code{K40/Ca44, err[K40/Ca44], Ca40/Ca44, err[Ca40/Ca44], rho}}
#' \item{\code{K40/Ca44, err[K40/Ca44], Ca40/Ca44, }
#'       \code{err[Ca40/Ca44], K40/Ca40, err[K40/Ca40]}}
#' }
#'
#' if \code{method='Rb-Sr'}, then \code{format} is one of either:
#'
#' \enumerate{
#' \item{\code{Rb87/Sr86, err[Rb87/Sr86], Sr87/Sr86, err[Sr87/Sr86] (, rho)}}
#' \item{\code{Rb, err[Rb], Sr, err[Sr], Sr87/Sr86, err[Sr87/Sr86]}}
#' }
#'
#' where \code{Rb} and \code{Sr} are in ppm
#'
#' if \code{method='Sm-Nd'}, then \code{format} is one of either:
#'
#' \enumerate{
#' \item{\code{Sm147/Nd144, err[Sm147/Nd144], Nd143/Nd144, err[Nd143/Nd144] (, rho)}}
#' \item{\code{Sm, err[Sm], Nd, err[Nd], Nd143/Nd144, err[Nd143/Nd144]}}
#' }
#'
#' where \code{Sm} and \code{Nd} are in ppm
#'
#' if \code{method='Re-Os'}, then \code{format} is one of either:
#'
#' \enumerate{
#' \item{\code{Re187/Os188, err[Re187/Os188], Os187/Os188, err[Os187/Os188] (, rho)}}
#' \item{\code{Re, err[Re], Os, err[Os], Os187/Os188, err[Os187/Os188]}}
#' }
#'
#' where \code{Re} and \code{Os} are in ppm
#'
#' if \code{method='Lu-Hf'}, then \code{format} is one of either:
#'
#' \enumerate{
#' \item{\code{Lu176/Hf177, err[Lu176/Hf177], Hf176/Hf177, err[Hf176/Hf177] (, rho)}}
#' \item{\code{Lu, err[Lu], Hf, err[Hf], Hf176/Hf177, err[Hf176/Hf177]}}
#' }
#'
#' where \code{Lu} and \code{Hf} are in ppm
#'
#' if \code{method='Th-U'}, then \code{format} is one of either:
#'
#' \enumerate{
#' \item{\code{X=8/2, err[X], Y=4/2, err[Y], Z=0/2, err[Z], rho[X,Y], rho[X,Z], rho[Y,Z]}}
#' \item{\code{X=2/8, err[X], Y=4/8, err[Y], Z=0/8, err[Z], rho[X,Y], rho[X,Z], rho[Y,Z]}}
#' \item{\code{X=8/2, err[X], Y=0/2, err[Y], rho[X,Y]}}
#' \item{\code{X=2/8, err[X], Y=0/8, err[Y], rho[X,Y]}}
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
#' or more U/Ca-ratios or U-concentration measurements (in ppm) and
#' their analytical uncertainties.}  }
#'
#' if \code{method='other'}, \code{x} is read as a table, unless
#' \code{format} is one of either:
#'
#' \describe{
#' \item{\code{radial} or \code{'average'}:}{\code{X, err[X]}}
#' \item{\code{regression}:}{\code{X, err[X], Y, err[Y], rho} \cr
#'       OR \code{X/Z, err[X/Z], Y/Z, err[Y/Z], X/Y, err[X/Y]}}
#' \item{\code{spectrum}:}{\code{f, X, err[X]}} 
#' }
#' 
#' @param ierr indicates whether the analytical uncertainties are
#'     reported as: \enumerate{
#' 
#' \item{1\eqn{\sigma} absolute uncertainties.}
#' \item{2\eqn{\sigma} absolute uncertainties.}
#' \item{1\eqn{\sigma} relative uncertainties (\%).}
#' \item{2\eqn{\sigma} relative uncertainties (\%).}
#' 
#' }
#'
#' @param U48 the initial \eqn{^{234}}U/\eqn{^{238}}U-activity ratio
#' @param Th0U8 the initial \eqn{^{230}}Th/\eqn{^{238}}U-activity ratio
#' @param Ra6U8 the initial \eqn{^{226}}Ra/\eqn{^{238}}U-activity ratio
#' @param Pa1U5 the initial \eqn{^{231}}Pa/\eqn{^{235}}U-activity ratio
#' 
#' @param Th02 2-element vector with the assumed initial
#'     \eqn{^{230}}Th/\eqn{^{232}}Th-ratio of the detritus and its
#'     standard error. Only used if \code{isochron==FALSE} and
#'     \code{detritus==2} 
#' @param Th02U48 9-element vector with the measured composition of
#'     the detritus, containing \code{X=0/8}, \code{sX}, \code{Y=2/8},
#'     \code{sY}, \code{Z=4/8}, \code{sZ}, \code{rXY}, \code{rXZ},
#'     \code{rYZ}.
#' 
#' @param ... optional arguments to the \code{read.csv} function
#' 
#' @seealso \code{\link{examples}}, \code{\link{settings}}
#' 
#' @return an object of class \code{UPb}, \code{PbPb}, \code{ArAr},
#'     \code{KCa}, \code{UThHe}, \code{ReOs}, \code{SmNd},
#'     \code{RbSr}, \code{LuHf}, \code{detritals},
#'     \code{fissiontracks}, \code{ThU} or \code{other}
#' @examples
#'
#' f1 <- system.file("UPb1.csv",package="IsoplotR")
#' file.show(f1) # inspect the contents of 'UPb1.csv'
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
#' f7 <- system.file("DZ.csv",package="IsoplotR")
#' d7 <- read.data(f7,method="detritals")
#' kde(d7)
#'
#' #  four 'other' files (LudwigMixture.csv, LudwigSpectrum.csv,
#' #  LudwigMean.csv, LudwigKDE.csv)
#' f8 <- system.file("LudwigMixture.csv",package="IsoplotR")
#' d8 <- read.data(f8,method="other")
#' radialplot(d8)
#'
#' @rdname read.data
#' @export
read.data <- function(x,...){ UseMethod("read.data",x) }
#' @rdname read.data
#' @export
read.data.default <- function(x,method='U-Pb',format=1,ierr=1,
                              U48=1,Th0U8=1,Ra6U8=1,Pa1U5=1,
                              Th02=c(0,0),Th02U48=c(0,0,1e6,0,0,0,0,0,0),...){
    X <- as.matrix(utils::read.table(x,sep=',',...))
    read.data.matrix(X,method=method,format=format,ierr=ierr,
                     U48=U48,Th0U8=Th0U8,Ra6U8=Ra6U8,Pa1U5=Pa1U5,
                     Th02=Th02,Th02U48=Th02U48)}
#' @rdname read.data
#' @export
read.data.data.frame <- function(x,method='U-Pb',format=1,ierr=1,
                                 U48=1,Th0U8=1,Ra6U8=1,Pa1U5=1,
                                 Th02=c(0,0),Th02U48=c(0,0,1e6,0,0,0,0,0,0),...){
    read.data.matrix(as.matrix(x),method=method,format=format,ierr=ierr,
                     U48=U48,Th0U8=Th0U8,Ra6U8=Ra6U8,Pa1U5=Pa1U5,
                     Th02=Th02,Th02U48=Th02U48,...)
}
#' @rdname read.data
#' @export
read.data.matrix <- function(x,method='U-Pb',format=1,ierr=1,
                             U48=1,Th0U8=1,Ra6U8=1,Pa1U5=1,
                             Th02=c(0,0),Th02U48=c(0,0,1e6,0,0,0,0,0,0),...){
    if (identical(method,'U-Pb')){
        out <- as.UPb(x,format=format,ierr=ierr,
                      U48=U48,Th0U8=Th0U8,Ra6U8=Ra6U8,Pa1U5=Pa1U5)
    } else if (identical(method,'Pb-Pb')){
        out <- as.PbPb(x,format=format,ierr=ierr)
    } else if (identical(method,'Ar-Ar')){
        out <- as.ArAr(x,format=format,ierr=ierr)
    } else if (identical(method,'K-Ca')){
        out <- as.KCa(x,format=format,ierr=ierr)
    } else if (identical(method,'Re-Os')){
        out <- as.ReOs(x,format=format,ierr=ierr)
    } else if (identical(method,'Rb-Sr')){
        out <- as.RbSr(x,format=format,ierr=ierr)
    } else if (identical(method,'Sm-Nd')){
        out <- as.SmNd(x,format=format,ierr=ierr)
    } else if (identical(method,'Lu-Hf')){
        out <- as.LuHf(x,format=format,ierr=ierr)
    } else if (identical(method,'Th-U')){
        out <- as.ThU(x,format=format,ierr=ierr,Th02=Th02,Th02U48=Th02U48)
    } else if (identical(method,'U-Th-He')){
        out <- as.UThHe(x,ierr=ierr)
    } else if (identical(method,'fissiontracks')){
        out <- as.fissiontracks(x,format=format,ierr=ierr)
    } else if (identical(method,'detritals')){
        out <- as.detritals(x)
    } else if (identical(method,'other')){
        out <- as.other(x,format=format,ierr=ierr)
    }
    out
}

as.UPb <- function(x,format=3,ierr=1,U48=1,Th0U8=1,Ra6U8=1,Pa1U5=1){
    out <- list()
    class(out) <- "UPb"
    out$x <- NA
    out$format <- format
    out$d <- diseq(U48=U48,Th0U8=Th0U8,Ra6U8=Ra6U8,Pa1U5=Pa1U5)
    nc <- ncol(x)
    nr <- nrow(x)
    if (is.numeric(x)) X <- x
    else X <- shiny2matrix(x,2,nr,nc)
    X <- errconvert(X,gc='U-Pb',format=format,ierr=ierr)
    cnames <- NULL
    if (format==1 & nc>4){
        cnames <- c('Pb207U235','errPb207U235',
                    'Pb206U238','errPb206U238','rhoXY')
    } else if (format==2 & nc>3){
        cnames <- c('U238Pb206','errU238Pb206',
                    'Pb207Pb206','errPb207Pb206','rhoXY')
        X <- read.XsXYsYrXY(X)
    } else if (format==3 & nc>5){
        cnames <- c('Pb207U235','errPb207U235',
                    'Pb206U238','errPb206U238',
                    'Pb207Pb206','errPb207Pb206',
                    'rhoXY','rhoYZ')
        if (nc > 7){
            rhoXY <- X[,7]
            rhoYZ <- X[,8]
            i <- which(is.na(rhoXY))
            j <- which(is.na(rhoYZ))
        } else if (nc == 7){
            rhoXY <- X[,7]
            i <- which(is.na(rhoXY))
            j <- 1:nrow(X)
            X <- cbind(X,0)
        } else {
            i <- 1:nrow(X)
            j <- 1:nrow(X)
            X <- cbind(X,0,0)
        }
        X[i,7] <- get.cor.75.68(X[i,1],X[i,2],X[i,3],X[i,4],X[i,5],X[i,6])
        X[j,8] <- get.cor.68.76(X[j,1],X[j,2],X[j,3],X[j,4],X[j,5],X[j,6])
    } else if (format==4 & nc>8){
        cnames <- c('Pb207U235','errPb207U235',
                    'Pb206U238','errPb206U238',
                    'Pb204U238','errPb204U238',
                    'rhoXY','rhoXZ','rhoYZ')
    } else if (format==5 & nc>8){
        cnames <- c('U238Pb206','errU238Pb206',
                    'Pb207Pb206','errPb207Pb206',
                    'Pb204Pb206','errPb204Pb206',
                    'rhoXY','rhoXZ','rhoYZ')
    } else if (format==6 & nc>11){
        cnames <- c('Pb207U235','errPb207U235',
                    'Pb206U238','errPb206U238',
                    'Pb204U238','errPb204U238',
                    'Pb207Pb206','errPb207Pb206',
                    'Pb204Pb207','errPb204Pb207',
                    'Pb204Pb206','errPb204Pb206')
    }
    out$x <- subset(X,select=1:length(cnames))
    colnames(out$x) <- cnames
    out
}
get.cor.75.68 <- function(Pb207U235,errPb207U235,
                          Pb206U238,errPb206U238,
                          Pb207Pb206,errPb207Pb206){
    get.cor.div(Pb207U235,errPb207U235,
                Pb206U238,errPb206U238,
                Pb207Pb206,errPb207Pb206)
}
get.cor.68.76 <- function(Pb207U235,errPb207U235,
                          Pb206U238,errPb206U238,
                          Pb207Pb206,errPb207Pb206){
    get.cor.mult(Pb206U238,errPb206U238,
                 Pb207Pb206,errPb207Pb206,
                 Pb207U235,errPb207U235)
}
get.cov.75.68 <- function(Pb207U235,errPb207U235,
                          Pb206U238,errPb206U238,
                          Pb207Pb206,errPb207Pb206){
    get.cov.div(Pb207U235,errPb207U235,
                Pb206U238,errPb206U238,
                Pb207Pb206,errPb207Pb206)
}
get.cov.75.48 <- function(Pb207U235,errPb207U235,
                          Pb204U238,errPb204U238,
                          Pb204Pb207,errPb204Pb207){
    get.cov.div(Pb207U235,errPb207U235,
                Pb204U238,errPb204U238,
                Pb204Pb207,errPb204Pb207)
}
get.cov.68.48 <- function(Pb206U238,errPb206U238,
                          Pb204U238,errPb204U238,
                          Pb204Pb206,errPb204Pb206){
    get.cov.div(Pb206U238,errPb206U238,
                Pb204U238,errPb204U238,
                Pb204Pb206,errPb204Pb206)
}
get.cov.76.86 <- function(Pb207Pb206,errPb207Pb206,
                          U238Pb206,errU238Pb206,
                          Pb207U235,errPb207U235){
    get.cov.div(Pb207Pb206,errPb207Pb206,
                U238Pb206,errU238Pb206,
                Pb207U235,errPb207U235)
}
get.cov.46.86 <- function(Pb204Pb206,errPb204Pb206,
                          U238Pb206,errU238Pb206,
                          Pb204U238,errPb204U238){
    get.cov.div(Pb204Pb206,errPb204Pb206,
                U238Pb206,errU238Pb206,
                Pb204U238,errPb204U238)
}
get.cov.46.68 <- function(Pb204Pb206,errPb204Pb206,
                          Pb206U238,errPb206U238,
                          Pb204U238,errPb204U238){
    get.cov.mult(Pb204Pb206,errPb204Pb206,
                 Pb206U238,errPb206U238,
                 Pb204U238,errPb204U238)
}
get.cov.46.76 <- function(Pb204Pb206,errPb204Pb206,
                          Pb207Pb206,errPb207Pb206,
                          Pb204Pb207,errPb204Pb207){
    get.cov.div(Pb204Pb206,errPb204Pb206,
                Pb207Pb206,errPb207Pb206,
                Pb204Pb207,errPb204Pb207)
}
get.cov.47.75 <- function(Pb204Pb207,errPb204Pb207,
                          Pb207U235,errPb207U235,
                          Pb204U238,errPb204U238){
    get.cov.mult(Pb204Pb207,errPb204Pb207,
                 Pb207U235,errPb207U235,
                 Pb204U238,errPb204U238)
}
as.PbPb <- function(x,format=1,ierr=1){
    out <- list()
    class(out) <- "PbPb"
    out$x <- NA
    out$format <- format
    nc <- ncol(x)
    nr <- nrow(x)
    if (is.numeric(x)) X <- x
    else X <- shiny2matrix(x,2,nr,nc)
    X <- errconvert(X,gc='Pb-Pb',format=format,ierr=ierr)
    cnames <- NULL
    if (format==1 & nc>4){
        cnames <- c('Pb206Pb204','errPb206Pb204',
                    'Pb207Pb204','errPb207Pb204','rho')
        out$x <- read.XsXYsYrXY(X)
    } else if (format==2 & nc>4) {
        cnames <- c('Pb204Pb206','errPb204Pb206',
                    'Pb207Pb206','errPb207Pb206','rho')
        out$x <- read.XsXYsYrXY(X)
    } else if (format==3 & nc>5){
        cnames <- c('Pb206Pb204','errPb206Pb204',
                    'Pb207Pb204','errPb207Pb204',
                    'Pb207Pb206','errPb207Pb206')
        out$x <- subset(X,select=1:length(cnames))
    }
    colnames(out$x) <- cnames
    out
}
as.ArAr <- function(x,format=3,ierr=1){
    out <- list()
    class(out) <- "ArAr"
    out$format <- format
    out$J <- errAdjust(x=as.numeric(x[2,1:2]),ierr=ierr)
    nc <- ncol(x)
    nr <- nrow(x)
    bi <- 4 # begin index
    X <- shiny2matrix(x,bi,nr,nc)
    X <- errconvert(X,gc='Ar-Ar',format=format,ierr=ierr)
    if (format==3 & nc>6){
        if (nc==7){
            out$x <- subset(X,select=1:7)
        } else {
            ns <- nr-bi+1 # number of samples
            out$x <- cbind(subset(X,select=1:6),1/ns)
        }
        colnames(out$x) <- c('Ar39Ar40','errAr39Ar40',
                             'Ar36Ar40','errAr36Ar40',
                             'Ar39Ar36','errAr39Ar36','Ar39')
    } else if (nc>3){
        if (nc>5){
            out$x <- subset(X,select=1:6)
        }
        if (nc==4) {
            X <- cbind(subset(X,select=1:4),0)
        }
        if (nc %in% c(4,5)){
            ns <- nr-bi+1 # number of samples
            out$x <- cbind(subset(X,select=1:5),1/ns)
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
as.KCa <- function(x,format=1,ierr=1){
    out <- list()
    class(out) <- "KCa"
    out$format <- format
    nc <- ncol(x)
    nr <- nrow(x)
    bi <- 2 # begin index
    X <- shiny2matrix(x,bi,nr,nc)
    X <- errconvert(X,gc='K-Ca',format=format,ierr=ierr)
    cnames <- c('K40Ca44','errK40Ca44','Ca40Ca44','errCa40Ca44','rho')
    if (format==1 & nc==4){
        out$x <- cbind(X,0)
    } else if (format==1 & nc>4){
        out$x <- subset(X,select=1:5)
    } else if (format==2 & nc>5){
        out$x <- subset(X,select=1:6)
        cnames <- c('K40Ca44','errK40Ca44',
                    'Ca40Ca44','errCa40Ca44',
                    'K40Ca40','errK40Ca40')
    }
    colnames(out$x) <- cnames
    out
}
as.RbSr <- function(x,format=2,ierr=1){
    cnames1 <- c('Rb87Sr86','errRb87Sr86',
                 'Sr87Sr86','errSr87Sr86','rho')
    cnames2 <- c('Rbppm','errRbppm','Srppm','errSrppm',
                 'Sr87Sr86','errSr87Sr86')
    as.PD(x,"RbSr",cnames1,cnames2,format,ierr)
}
as.ReOs <- function(x,format=2,ierr=1){
    cnames1 <- c('Re187Os188','errRe187Os188',
                 'Os187Os188','errOs187Os188','rho')
    cnames2 <- c('Reppm','errReppm','Osppm','errOsppm',
                 'Os187Os188','errOs187Os188')
    as.PD(x,"ReOs",cnames1,cnames2,format,ierr)
}
as.SmNd <- function(x,format=2,ierr=1){
    cnames1 <- c('Sm143Nd144','errSm143Nd144',
                 'Nd143Nd144','errNd143Nd144','rho')
    cnames2 <- c('Smppm','errSmppm','Ndppm','errNdppm',
                 'Nd143Nd144','errNd143Nd144')
    as.PD(x,"SmNd",cnames1,cnames2,format,ierr)
}
as.LuHf <- function(x,format=2,ierr=1){
    cnames1 <- c('Lu176Hf177','errLu176Hf177',
                 'Hf176Hf177','errHf176Hf177','rho')
    cnames2 <- c('Luppm','errLuppm','Hfppm','errHfppm',
                 'Hf176Hf177','errHf176Hf177')
    as.PD(x,"LuHf",cnames1,cnames2,format,ierr)
}
as.PD <- function(x,classname,colnames1,colnames2,format,ierr){
    out <- list()
    class(out) <- c(classname,'PD')
    out$x <- NA
    out$format <- format
    nc <- ncol(x)
    nr <- nrow(x)
    if (is.numeric(x)) X <- x
    else X <- shiny2matrix(x,2,nr,nc)
    X <- errconvert(X,gc='PD',format=format,ierr=ierr)
    if (format==1 & nc>3){
        out$x <- read.XsXYsYrXY(X)
        colnames(out$x) <- colnames1
    } else if (format==2 & nc>5){
        out$x <- subset(X,select=1:6)
        colnames(out$x) <- colnames2
    }
    out
}
as.ThU <- function(x,format=1,ierr=1,Th02=c(0,0),Th02U48=c(0,0,1e6,0,0,0,0,0,0)){
    out <- list()
    class(out) <- "ThU"
    out$x <- NA
    out$format <- format
    out$Th02 <- Th02
    out$Th02U48 <- Th02U48
    nc <- ncol(x)
    nr <- nrow(x)
    if (is.numeric(x)) X <- x
    else X <- shiny2matrix(x,2,nr,nc)
    X <- errconvert(X,gc='Th-U',format=format,ierr=ierr)
    cnames <- NULL
    if (format==1 & nc>8){
        cnames <- c('U238Th232','errU238Th232',
                    'U234Th232','errU234Th232',
                    'Th230Th232','errTh230Th232',
                    'rhoXY','rhoXZ','rhoYZ')
        out$x <- subset(X,select=1:9)
    } else if (format==2 & nc>8) {
        cnames <- c('Th232U238','errTh232U238',
                    'U234U238','errU234U238',
                    'Th230U238','errTh230U238',
                    'rhoXY','rhoXZ','rhoYZ')
        out$x <- subset(X,select=1:9)
    } else if (format==3 & nc>3) {
        if (nc==4) X <- cbind(subset(X,select=1:4),0)
        cnames <- c('U238Th232','errU238Th232',
                    'Th230Th232','errTh230Th232',
                    'rho')
        out$x <- subset(X,select=1:5)
    } else if (format==4 & nc>3) {
        if (nc==4) X <- cbind(subset(X,select=1:4),0)
        cnames <- c('Th232U238','errTh232U238',
                    'Th230U238','errTh230U238',
                    'rho')
        out$x <- subset(X,select=1:5)
    }
    out$x <- subset(X,select=1:length(cnames))
    colnames(out$x) <- cnames
    out
}
as.UThHe <- function(x,ierr=1){
    nc <- ncol(x)
    nr <- nrow(x)
    if (is.numeric(x)) X <- x
    else X <- shiny2matrix(x,2,nr,nc)
    X[X<=0] <- NA
    X <- errconvert(X,gc='U-Th-He',ierr=ierr)
    if (nc>5) cnames <- c('He','errHe','U','errU','Th','errTh')
    if (nc>7) cnames <- c(cnames,'Sm','errSm')
    out <- subset(X,select=1:length(cnames))
    colnames(out) <- cnames
    class(out) <- append("UThHe",class(out))
    out
}
as.fissiontracks <- function(x,format=1,ierr=1){
    nr <- nrow(x)
    nc <- ncol(x)
    out <- list()
    class(out) <- "fissiontracks"
    out$format <- format
    if (format==1){
        out$zeta <- errAdjust(as.numeric(x[2,1:2]),ierr=ierr)
        out$rhoD <- errAdjust(as.numeric(x[4,1:2]),ierr=ierr)
        X <- shiny2matrix(x,6,nr,nc)
        out$x <- subset(X,select=1:2)
        colnames(out$x) <- c('Ns','Ni')
    } else {
        if (format==2){
            out$zeta <- errAdjust(as.numeric(x[2,1:2]),ierr=ierr)
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
        errcols <- seq(from=4,to=nc,by=2)-2
        for (i in 1:ns){
            UsU <- as.numeric(x[i+si-1,3:nc])
            UsU <- errAdjust(UsU,i=errcols,ierr=ierr)
            out$U[[i]] <- UsU[1]
            out$sU[[i]] <- UsU[2]
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
    if (is.numeric(x)) X <- x
    else X <- shiny2matrix(x,2,nr,nc)
    colnames(X) <- snames
    for (sname in snames){
        out[[sname]] = X[!is.na(X[,sname]),sname]
    }
    out
}
as.other <- function(x,format='generic',ierr=1){
    nc <- ncol(x)
    has.header <- is.na(suppressWarnings(as.numeric(x[1,1])))
    if (has.header) x <- x[-1,]
    X <- matrix(as.numeric(x),ncol=nc)
    errconvert(X,gc='other',format=format,ierr=ierr)
}

# x = a numerical vector, br = length of the preamble with parameters
# nr = number of rows, nc = number of columns
shiny2matrix <- function(x,br,nr,nc){
    suppressWarnings(
        return(matrix(as.numeric(x[(br:nr),]),nr-br+1,nc))
    )
}

# for data of class UPb, PbPb, PD (including LuHf, SmNd, RbSr, and ReOs) 
read.XsXYsYrXY <- function(x){
    nc <- ncol(x)
    if (nc == 4){
        out <- cbind(x,0)
    } else {
        i <- which(is.na(x[,5]))
        x[i,5] <- 0
        out <- subset(x,select=1:5)
    }
    out
}

errconvert <- function(x,gc='U-Pb',format=1,ierr=1){
    if (ierr==1){ return(x) }
    else { out <- x }
    i <- getErrCols(gc,format,ierr)
    errAdjust(x,i,ierr)
}

errAdjust <- function(x,i=2,ierr=1){
    out <- x
    if (is.vector(x)){
        if (ierr==2){
            out[i] <- x[i]/2
        } else if (ierr==3){
            out[i] <- x[i]*x[i-1]/100
        } else if (ierr==4){
            out[i] <- x[i]*x[i-1]/200
        }
    } else {
        if (ierr==2){
            out[,i] <- x[,i]/2
        } else if (ierr==3){
            out[,i] <- x[,i]*x[,i-1]/100
        } else if (ierr==4){
            out[,i] <- x[,i]*x[,i-1]/200
        }
    }
    out
}

getErrCols <- function(gc,format=NA,ierr=1){
    UPb12 = (gc=='U-Pb' && format%in%(1:2))
    UPb345 = (gc=='U-Pb' && format%in%(3:5))
    UPb6 = (gc=='U-Pb' && format==6)
    PbPb12 = (gc=='Pb-Pb' && format%in%(1:2))
    PbPb3 = (gc=='Pb-Pb' && format==3)
    ArAr12 = (gc=='Ar-Ar' && format%in%(1:2))
    ArAr3 = (gc=='Ar-Ar' && format==3)
    KCa1 = (gc=='K-Ca' && format==1)
    KCa2 = (gc=='K-Ca' && format==2)
    PD1 = (gc=='PD' && format==1)
    PD2 = (gc=='PD' && format==2)
    UThHe = (gc=='U-Th-He')
    ThU12 = (gc=='Th-U' && format<3)
    ThU34 = (gc=='Th-U' && format>2)
    radial = (gc=='other' && identical(format,'radial'))
    regression = (gc=='other' && identical(format,'regression'))
    spectrum = (gc=='other' && identical(format,'spectrum'))
    average = (gc=='other' && identical(format,'average'))
    if (UPb12 | PbPb12 | ArAr12 | KCa1 | PD1 | ThU34 | regression){
        cols = c(2,4)
    } else if (UPb345 | PbPb3 | ArAr3 | KCa2 | PD2 | UThHe | ThU12){
        cols = c(2,4,6)
    } else if (UPb6){
        cols = seq(from=2,to=12,by=2)
    } else if (radial || average){
        cols = 2
    } else if (spectrum){
        cols = 3
    } else {
        cols = NULL
    }
    cols
}
