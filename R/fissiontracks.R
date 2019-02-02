fissiontrack.age <- function(x,i=NA,sigdig=NA,exterr=TRUE,mineral=NA){
    if (x$format < 2){
        out <- EDM.age(x,i,sigdig,exterr)
    } else if (x$format > 1){
        if (x$format == 3) {
            if (is.na(mineral)) mineral <- 'apatite'
            x$zeta <- get.absolute.zeta(mineral);
        }
        out <- ICP.age(x,i,sigdig,exterr)
    }
    out
}

get.absolute.zeta <- function(mineral){
    R <- iratio('U238U235')[1]
    MM <- imass('U')[1]
    qap <- etchfact(mineral)
    L <- tracklength(mineral)
    Lf <- lambda('fission')[1]
    dens <- mindens(mineral)
    Na <- 6.02214e23
    zeta <- 4*(1+R)*MM*1e18/(Na*Lf*qap*L*dens*R)
    c(zeta,0)
}

#' Calculate the zeta calibration coefficient for fission track dating
#'
#' Determines the zeta calibration constant of a fission track dataset
#' (EDM or LA-ICP-MS) given its true age and analytical uncertainty.
#'
#' @details
#' The fundamental fission track age is given by:
#'
#' \eqn{t = \frac{1}{\lambda_{238}}
#' \ln\left(1 + \frac{\lambda_{238}}{\lambda_f} \frac{2 N_s}{[^{238}U]A_sL}\right)
#' } (eq.1)
#'
#' where \eqn{N_s} is the number of spontaneous fission tracks
#' measured over an area \eqn{A_s}, \eqn{[^{238}U]} is the
#' \eqn{^{238}}U-concentration in atoms per unit volume,
#' \eqn{\lambda_f} is the fission decay constant, \eqn{L} is the
#' etchable fission track length, and the factor 2 is a geometric
#' factor accounting for the fact that etching reveals tracks from
#' both above and below the internal crystal surface.  Two analytical
#' approaches are used to measure \eqn{[^{238}U]}: neutron activation
#' and LAICPMS. The first approach estimates the
#' \eqn{^{238}}U-concentration indirectly, using the induced fission
#' of neutron-irradiated \eqn{^{235}}U as a proxy for the
#' \eqn{^{238}}U. In the most common implementation of this approach,
#' the induced fission tracks are recorded by an external detector
#' made of mica or plastic that is attached to the polished grain
#' surface (Fleischer and Hart, 1972; Hurford and Green, 1983). The
#' fission track age equation then becomes:
#'
#' \eqn{t = \frac{1}{\lambda_{238}}
#' \ln\left(1 + \frac{\lambda_{238}\zeta\rho_d}{2}\frac{N_s}{N_i}\right)
#' } (eq.2)
#'
#' where \eqn{N_i} is the number of induced fission tracks counted in
#' the external detector over the same area as the spontaneous tracks,
#' \eqn{\zeta} is a `zeta'-calibration factor that incorporates both
#' the fission decay constant and the etchable fission track length,
#' and \eqn{\rho_d} is the number of induced fission tracks per unit
#' area counted in a co-irradiated glass of known
#' U-concentration. \eqn{\rho_d} allows the \eqn{\zeta}-factor to be
#' `recycled' between irradiations.
#'
#' LAICPMS is an alternative means of determining the
#' \eqn{^{238}}U-content of fission track samples without the need for
#' neutron irradiation. The resulting U-concentrations can be plugged
#' directly into the fundamental age equation (eq.1). but this is
#' limited by the accuracy of the U-concentration measurements, the
#' fission track decay constant and the etching and counting
#' efficiencies. Alternatively, these sources of bias may be removed
#' by normalising to a standard of known fission track age and
#' defining a new `zeta' calibration constant \eqn{\zeta_{icp}}:
#'
#' \eqn{
#' t = \frac{1}{\lambda_{238}}
#' \ln\left( 1 + \frac{\lambda_{238}\zeta_{icp}}{2} \frac{N_s}{[{}^{238}U] A_s} \right)
#' }(eq.3)
#'
#' where \eqn{[{}^{238}U]} may either stand for the
#' \eqn{^{238}}U-concentration (in ppm) \emph{or} for the U/Ca (for
#' apatite) or U/Si (for zircon) ratio measurement (Vermeesch, 2017).
#'
#' @param x an object of class \code{fissiontracks}
#' @param tst a two-element vector with the true age and its standard
#'     error
#' @param exterr logical flag indicating whether the external
#'     uncertainties associated with the age standard or the dosimeter
#'     glass (for the EDM) should be accounted for when propagating
#'     the uncertainty of the zeta calibration constant.
#' @param update logical flag indicating whether the function should
#'     return an updated version of the input data, or simply return a
#'     two-element vector with the calibration constant and its
#'     standard error.
#' @param sigdig number of significant digits
#' @return an object of class \code{fissiontracks} with an updated
#'     \code{x$zeta} value
#' @seealso \code{\link{age}}
#' @examples
#' data(examples)
#' print(examples$FT1$zeta)
#' FT <- set.zeta(examples$FT1,tst=c(250,5))
#' print(FT$zeta)
#'
#' @references
#' Fleischer, R. and Hart, H. Fission track dating: techniques and
#' problems. In Bishop, W., Miller, J., and Cole, S., editors,
#' Calibration of Hominoid Evolution, pages 135-170. Scottish Academic
#' Press Edinburgh, 1972.
#'
#' Hurford, A. J. and Green, P. F. The zeta age calibration of
#' fission-track dating. Chemical Geology, 41:285-317, 1983.
#'
#' Vermeesch, P., 2017. Statistics for LA-ICP-MS based fission track
#' dating. Chemical Geology, 456, pp.19-27.
#' @export
set.zeta <- function(x,tst,exterr=TRUE,update=TRUE,sigdig=2){
    N <- length(x$Ns)
    L8 <- lambda('U238')[1]
    tt <- tst[1]
    if (exterr) st <- tst[2]
    else st <- 0
    if (x$format==1){
        Ns <- sum(x$x[,'Ns'])
        Ni <- sum(x$x[,'Ni'])
        rhoD <- x$rhoD
        if (!exterr) rhoD[2] <- 0
        zeta <- 2e6*(exp(L8*tt)-1)/(L8*rhoD[1]*Ns/Ni)
        zetaErr <- zeta * sqrt( (L8*exp(L8*tt)*st/(exp(L8*tt)-1))^2 +
                                (rhoD[2]/rhoD[1])^2 + 1/Ns + 1/Ni )
    } else {
        Ns <- sum(x$Ns)
        UsU <- get.UsU(x)
        UA <- sum(UsU[,1]*x$A)
        UAerr <- sqrt( sum(UsU[,2]*x$A)^2 )
        zeta <- 2*UA*(exp(L8*tt)-1)/(L8*Ns)
        zetaErr <- zeta * sqrt( ((L8*exp(L8*tt)*st)/(exp(L8*tt)-1))^2 +
                                1/Ns + (UAerr/UA)^2 )
    }
    zsz <- roundit(zeta,zetaErr,sigdig=sigdig)
    if (update){
        out <- x
        out$zeta <- zsz
    } else {
        out <- matrix(zsz,1,2)
        colnames(out) <- c('zeta','s[zeta]')
    }
    out
}

ICP.age <- function(x,i=NA,sigdig=NA,exterr=TRUE){
    ngrains <- length(x$Ns)
    tt <- rep(NA,ngrains)
    st <- rep(NA,ngrains)
    if (exterr){
        zeta <- x$zeta
    } else {
        zeta <- c(x$zeta[1],0)
    }
    ipos <- which(x$Ns>0)
    izero <- which(x$Ns<1)
    UsU <- get.UsU(x)
    # first calculate the ages of the non-zero track data:
    for (i in ipos){
        tst <- get.ICP.age(x$Ns[i],x$A[i],UsU[i,],zeta)
        tt[i] <- tst[1]
        st[i] <- tst[2]
    }
    # then use the equivalent induced track approach:
    for (i in izero){
        rho <- UsU[i,'U']/(x$A[i]*UsU[i,'sU']^2)
        Ni <- x$A[i]*UsU[i,'U']*rho
        tst <- get.EDM.age(x$Ns[i],Ni,zeta*1e6,c(rho,0))
        tt[i] <- tst[1]
        st[i] <- tst[2]
    }
    out <- roundit(tt,st,sigdig=sigdig)
    colnames(out) <- c('t','s[t]')
    out
}

get.UsU <- function(x){
    Aicp <- pi*(x$spotSize/2)^2
    n <- length(x$U)
    nspots <- length(stats::na.omit(unlist(x$U)))
    do.average <- (nspots>n)
    out <- matrix(0,n,2)
    colnames(out) <- c('U','sU')
    m <- rep(0,n)
    if (do.average) {
        uhat <- rep(0,n)
        num <- 0
        den <- 0
    }
    for (j in 1:n){
        if (do.average){
            Uj <- stats::na.omit(x$U[[j]])
            m[j] <- length(Uj) # spots per grain
            uhat[j] <- mean(log(Uj))
            num <- num + sum((log(Uj)-uhat[j])^2)
            den <- den + m[j] - 1
        } else {
            out[j,'U'] <- stats::na.omit(x$U[[j]])
            out[j,'sU'] <- stats::na.omit(x$sU[[j]])
        }
    }
    if (do.average){
        out[,'U'] <- exp(uhat)
        vhat <- rep(num/den,n)
        for (j in 1:n){
            suhat <- x$sU[[j]]/x$U[[j]]
            vhat[j] <- vhat[j]*(1-m[j]*Aicp/x$A[j])^2 +
                       sum(suhat^2,na.rm=TRUE)*(Aicp/x$A[j])^2
        }
        out[,'sU'] <- exp(uhat)*sqrt(vhat)
    }
    out
}

EDM.age <- function(x,i=NA,sigdig=2,exterr=TRUE){
    ns <- nrow(x$x)
    out <- matrix(0,ns,2)
    colnames(out) <- c('t','s[t]')
    if (exterr){
        zeta <- x$zeta
        rhoD <- x$rhoD
    } else {
        zeta <- c(x$zeta[1],0)
        rhoD <- c(x$rhoD[1],0)
    }
    for (j in 1:ns){
        tt <- get.EDM.age(x$x[j,'Ns'],x$x[j,'Ni'],zeta,rhoD)
        out[j,] <- roundit(tt[1],tt[2],sigdig=sigdig)
    }
    if (!is.na(i)) out <- out[i,]
    out
}

# zeta and rhoD are two-element vectors
get.EDM.age <- function(Ns,Ni,zeta,rhoD=c(1,0)){
    L8 <- lambda('U238')[1]
    if (Ns<1){
        Ns <- Ns+0.5
        Ni <- Ni+0.5
    }
    tt <- log(1+0.5*L8*(zeta[1]/1e6)*rhoD[1]*(Ns/Ni))/L8
    st <- tt*sqrt(1/Ns + 1/Ni + (rhoD[2]/rhoD[1])^2 + (zeta[2]/zeta[1])^2)
    c(tt,st)
}
# zeta is a two-element vector
get.ICP.age <- function(Ns,A,UsU,zeta){
    L8 <- lambda('U238')[1]
    tt <- log(1+L8*zeta[1]*Ns/(2*UsU[1]*A))/L8
    st <- tt * sqrt(1/Ns + (zeta[2]/zeta[1])^2 + (UsU[2]/UsU[1])^2)
    c(tt,st)
}
