#' @title Set up U-series disequilibrium correction for U-Pb
#'     geochronology
#' 
#' @description The U-Pb method conventionally assumes initial secular
#'     equilibrium of all the intermediate daughters of the
#'     \eqn{{}^{238}}U-\eqn{{}^{206}}Pb and
#'     \eqn{{}^{235}}U-\eqn{{}^{207}}Pb decay chains.  Violation of
#'     this assumption may produce inaccurate results.  \code{diseq}
#'     sets up initial disequilibrium parameters that are subsequently
#'     passed on to the \code{read.data} function for incorporation in
#'     other functions.
#'
#' @details
#' There are three ways to correct for the initial disequilibrium
#' between the activity of \eqn{{}^{238}}U, \eqn{{}^{234}}U,
#' \eqn{{}^{230}}Th, and \eqn{{}^{226}}Ra; or between \eqn{{}^{235}}U
#' and \eqn{{}^{231}}Pa:
#'
#' \enumerate{
#' 
#' \item{Specify the assumed initial activity ratios and calculate how
#' much excess \eqn{{}^{206}}Pb and \eqn{{}^{207}}Pb these would have
#' produced.}
#' 
#' \item{Measure the current activity ratios to infer the initial
#' ratios.  This approach only works for young samples.}
#' 
#' \item{The initial \eqn{{}^{230}}Th/\eqn{{}^{238}}U activity ratio
#' can also be estimated by providing the Th/U-ratio of the magma.}
#' 
#' }
#'
#' @param U48 a list containing seven items (\code{x}, \code{sx},
#'     \code{m}, \code{M}, \code{x0}, \code{sd} and \code{option})
#'     specifying the \eqn{{}^{234}}U/\eqn{{}^{238}}U disequilibrium.
#'
#' If \code{option=0}, then \code{x} and \code{sx} are ignored and no
#'     disequilibrium correction is applied.
#'
#' If \code{option=1}, then \code{x} contains the initial
#'     \eqn{{}^{234}}U/\eqn{{}^{238}}U ratio and \code{sx} its
#'     standard error.
#'
#' If \code{option=2}, then \code{x} contains the measured
#'     \eqn{{}^{234}}U/\eqn{{}^{238}}U ratio and \code{sx} its
#'     standard error.
#'
#' \code{m}, \code{M} specify the minimum and maximum search limits of
#' the initial \eqn{{}^{234}}U/\eqn{{}^{238}}U activity ratio.
#'
#' \code{x0} and \code{sd} specify the mean and standard deviation of
#' the prior distribution for the the initial
#' \eqn{{}^{234}}U/\eqn{{}^{238}}U activity ratio.
#'
#' @param ThU a list containing seven items (\code{x}, \code{sx},
#'     \code{m}, \code{M}, \code{x0}, \code{sd} and \code{option})
#'     specifying the \eqn{{}^{230}}Th/\eqn{{}^{238}}U disequilibrium.
#'
#' If \code{option=0}, then \code{x} and \code{sx} are ignored and no
#'     disequilibrium correction is applied.
#'
#' If \code{option=1}, then \code{x} contains the initial
#'     \eqn{{}^{230}}Th/\eqn{{}^{238}}U ratio and \code{sx} its
#'     standard error.
#'
#' If \code{option=2}, then \code{x} contains the measured
#'     \eqn{{}^{230}}Th/\eqn{{}^{238}}U ratio and \code{sx} its
#'     standard error.
#'
#' If \code{option=3}, then \code{x} contains the measured Th/U ratio
#'     of the magma (assumed or determined from the whole rock or
#'     volcanic glass) and \code{sx} its standard error. This only
#'     applies for Th-bearing U-Pb data formats 7 and 8.
#'
#' \code{m}, \code{M}, \code{x0} and \code{sd} are analogous to the
#' eponymous settings for \code{ThU}.
#'
#' @param RaU a list containing seven items (\code{x}, \code{sx},
#'     \code{m}, \code{M}, \code{x0}, \code{sd} and \code{option})
#'     specifying the \eqn{{}^{226}}Ra/\eqn{{}^{238}}U disequilibrium.
#'
#' If \code{option=0}, then \code{x} and \code{sx} are ignored and no
#'     disequilibrium correction is applied.
#'
#' If \code{option=1}, then \code{x} contains the initial
#'     \eqn{{}^{226}}Ra/\eqn{{}^{238}}U ratio and \code{sx} its
#'     standard error.
#' 
#' \code{m}, \code{M}, \code{x0} and \code{sd} are analogous to the
#' eponymous settings for \code{ThU}.
#' 
#' @param PaU a list containing seven items (\code{x}, \code{sx},
#'     \code{m}, \code{M}, \code{x0}, \code{sd} and \code{option})
#'     specifying the \eqn{{}^{231}}Pa/\eqn{{}^{235}}U disequilibrium.
#'
#' If \code{option=0}, then \code{x} and \code{sx} are ignored and no
#'     disequilibrium correction is applied.
#'
#' If \code{option=1}, then \code{x} contains the initial
#'     \eqn{{}^{231}}Pa/\eqn{{}^{235}}U ratio and \code{sx} its
#'     standard error.
#'
#' \code{m}, \code{M}, \code{x0} and \code{sd} are analogous to the
#' eponymous settings for \code{ThU}.
#' 
#' @param buffer small amount of padding to avoid singularities in the
#'     prior distribution, which uses a logistic transformation:
#' \eqn{y = \ln\left(\frac{x-m+buffer}{M+buffer-x}\right)}
#' 
#' @return a list with the following items:
#'
#' \describe{
#' 
#' \item{U48, ThU, RaU, PaU}{the same as the corresponding input arguments}
#' 
#' \item{equilibrium}{a boolean flag indicating whether
#' \code{option=TRUE} and/or \code{x=1} for all activity ratios}
#'
#' \item{Q}{the eigenvectors of the disequilibrium matrix exponential}
#'
#' \item{Qinv}{the inverse of \code{Q}}
#'
#' \item{L}{a named vector of all the relevant decay constants}
#' 
#' }
#'
#' @seealso \code{\link{mclean}}, \code{\link{concordia}},
#'     \code{\link{ludwig}}
#' @examples
#' d <- diseq(U48=list(x=0,option=1),ThU=list(x=2,option=1),
#'            RaU=list(x=2,option=1),PaU=list(x=2,option=1))
#' fn <- system.file("diseq.csv",package="IsoplotR")
#' UPb <- read.data(fn,method='U-Pb',format=2,d=d)
#' concordia(UPb,type=2,xlim=c(0,700),ylim=c(0.05,0.5))
#' @export
diseq <- function(U48=list(x=1,sx=0,option=0,m=0,M=20,x0=1,sd=10),
                  ThU=list(x=1,sx=0,option=0,m=0,M=20,x0=1,sd=10),
                  RaU=list(x=1,sx=0,option=0,m=0,M=20,x0=1,sd=10),
                  PaU=list(x=1,sx=0,option=0,m=0,M=20,x0=1,sd=10),
                  buffer=1e-5){
    patch <- function(aratio){
        out <- aratio
        if (is.null(aratio$x)) out$x <- 1
        if (is.null(aratio$sx)) out$sx <- 0
        if (is.null(aratio$option)) out$option <- 0
        if (is.null(aratio$m)) out$m <- 0
        if (is.null(aratio$M)) out$M <- 20
        if (is.null(aratio$x0)) out$x0 <- 1
        if (is.null(aratio$sd)) out$sd <- 10
        if ((out$x0<out$m)||(out$x0>out$M)) out$x0 <- (out$m+out$M)/2
        out
    }
    out <- list()
    class(out) <- 'diseq'
    out$U48 <- patch(U48)
    out$ThU <- patch(ThU)
    out$RaU <- patch(RaU)
    out$PaU <- patch(PaU)
    out$buffer <- buffer
    out$equilibrium <- check_equilibrium(d=out)
    l38 <- settings('lambda','U238')[1]
    l34 <- settings('lambda','U234')[1]*1000
    l30 <- settings('lambda','Th230')[1]*1000
    l26 <- settings('lambda','Ra226')[1]*1000
    l35 <- settings('lambda','U235')[1]
    l31 <- settings('lambda','Pa231')[1]*1000
    out$L <- c(l38,l34,l30,l26,0,l35,l31,0)
    names(out$L) <- c('U238','U234','Th230','Ra226','Pb206','U235','Pa231','Pb207')
    out$Q <- matrix(0,8,8)
    out$Q[1,1] <- -((l26-l38)*(l30-l38)*(l34-l38))/(l26*l30*l34)
    out$Q[2,1] <- -(l38*(l26-l38)*(l30-l38))/(l26*l30*l34)
    out$Q[2,2] <- -((l26-l34)*(l30-l34))/(l26*l30)
    out$Q[3,1] <- -(l38*(l26-l38))/(l26*l30)
    out$Q[3,2] <- -(l34*(l26-l34))/(l26*l30)
    out$Q[3,3] <- -(l26-l30)/l26
    out$Q[4,1] <- -l38/l26
    out$Q[4,2] <- -l34/l26
    out$Q[4,3] <- -l30/l26
    out$Q[4,4] <- -1
    out$Q[5,1:5] <- 1
    out$Q[6,6] <- -(l31-l35)/l31
    out$Q[7,6] <- -l35/l31
    out$Q[7,7] <- -1
    out$Q[8,6:8] <- 1
    out$Qinv <- matrix(0,8,8)
    out$Qinv[1,1] <- -(l26*l30*l34)/((l26-l38)*(l30-l38)*(l34-l38))
    out$Qinv[2,1] <- (l26*l30*l38)/((l26-l34)*(l30-l34)*(l34-l38))
    out$Qinv[2,2] <- -(l26*l30)/((l26-l34)*(l30-l34))
    out$Qinv[3,1] <- -(l26*l34*l38)/((l26-l30)*(l30-l34)*(l30-l38))
    out$Qinv[3,2] <- (l26*l34)/((l26-l30)*(l30-l34))
    out$Qinv[3,3] <- -l26/(l26-l30)
    out$Qinv[4,1] <- (l30*l34*l38)/((l26-l30)*(l26-l34)*(l26-l38))
    out$Qinv[4,2] <- -(l30*l34)/((l26-l30)*(l26-l34))
    out$Qinv[4,3] <- l30/(l26-l30)
    out$Qinv[4,4] <- -1
    out$Qinv[5,1:5] <- 1
    out$Qinv[6,6] <- -l31/(l31-l35)
    out$Qinv[7,6] <- l35 /(l31-l35)
    out$Qinv[7,7] <- -1
    out$Qinv[8,6:8] <- 1
    out
}

# infer initial activity ratios from Th/U measurements
copy_diseq <- function(x,d=diseq()){
    out <- d
    if (d$ThU$option==3){
        if (x$format>6){
            U <- settings('iratio','U238U235')[1]
            out$ThU$x <- (x$x[,'Th232U238']/d$ThU$x)*U/(1+U)
        }
        out$ThU$option <- 1
    }
    out
}

# to infer initial 38,34,30,35 from measured 34/38 and 30/38
mexp_8405 <- function(){
    l38 <- settings('lambda','U238')[1]
    l34 <- settings('lambda','U234')[1]*1000
    l30 <- settings('lambda','Th230')[1]*1000
    l35 <- settings('lambda','U235')[1]
    out <- list()
    out$L <- c(l38,l34,l30,l35)
    names(out$L) <- c('U238','U234','Th230','U235')
    out$Q <- diag(1,4,4)
    out$Q[1,1] <- (l30-l38)*(l34-l38)/(l34*l38)
    out$Q[2,1] <- (l30-l38)/l34
    out$Q[2,2] <- (l30-l34)/l34
    out$Q[3,1:2] <- 1
    out$Qinv <- diag(1,4,4)
    out$Qinv[1,1] <- (l34*l38)/((l30-l38)*(l34-l38))
    out$Qinv[2,1] <- -(l34*l38)/((l30-l34)*(l34-l38))
    out$Qinv[2,2] <- l34/(l30-l34)
    out$Qinv[3,1] <- (l34*l38)/((l30-l34)*(l30-l38))
    out$Qinv[3,2] <- -l34/(l30-l34)
    out
}

# to infer initial 38,34,35 from measured 34/38
mexp_845 <- function(nratios=3){
    l38 <- settings('lambda','U238')[1]
    l34 <- settings('lambda','U234')[1]*1000
    l30 <- settings('lambda','Th230')[1]*1000
    l26 <- settings('lambda','Ra226')[1]*1000
    l35 <- settings('lambda','U235')[1]
    out <- list()
    out$L <- c(l38,l34,l35)
    names(out$L) <- c('U238','U234','U235')
    out$Q <- diag(1,3,3)
    out$Q[1,1] <- (l34-l38)/l38
    out$Q[2,1] <- 1
    out$Qinv <- diag(1,3,3)
    out$Qinv[1,1] <- l38/(l34-l38)
    out$Qinv[2,1] <- -l38/(l34-l38)
    out
}

reverse <- function(tt,mexp,nt){
    n0 <- (mexp$Q %*% diag(exp(mexp$L*tt)) %*% mexp$Qinv %*% nt)
    rownames(n0) <- rownames(nt)
    n0
}
forward <- function(tt,d=diseq(),derivative=0){
    if (derivative==1){
        nt <- (d$Q %*% diag(exp(-d$L)) %*%
               diag(exp(-d$L*tt)) %*% d$Qinv %*% d$n0)
    } else if (derivative==2){
        nt <- (d$Q %*% diag(exp(-d$L)) %*% diag(exp(-d$L)) %*%
               diag(exp(-d$L*tt)) %*% d$Qinv %*% d$n0)
    } else {
        nt <- (d$Q %*% diag(exp(-d$L*tt)) %*% d$Qinv %*% d$n0)
    }
    rownames(nt) <- names(d$L)
    nt
}

check_equilibrium <- function(d=diseq()){
    checkratio <- function(r){
        r$option == 0 || all(r$x==1) & all(r$sx==0)
    }
    U48 <- checkratio(d$U48)
    ThU <- checkratio(d$ThU)
    RaU <- checkratio(d$RaU)
    PaU <- checkratio(d$PaU)
    U48 & ThU & RaU & PaU
}

measured_disequilibrium <- function(d=diseq()){
    measured <- (d$U48$option==2 | d$ThU$option==2)
    !d$equilibrium & measured
}

#' @title Predict disequilibrium concordia_compositions
#' 
#' @description
#' Returns the predicted \eqn{{}^{206}}Pb/\eqn{{}^{238}}U and
#' \eqn{{}^{207}}Pb/\eqn{{}^{235}}U ratios for any given time with or
#' without initial U-series disequilibrium.
#'
#' @details
#' U decays to Pb in 14 (for \eqn{{}^{238}}U) or 11/12 (for
#' \eqn{{}^{235}}U) steps. Conventional U-Pb geochronology assumes
#' that secular equilibrium between all the short lived intermediate
#' daughters was established at the time of isotopic closure. Under
#' this assumption, the relative abundances of those intermediate
#' daughters can be neglected and the age equation reduces to a simple
#' function of the measured Pb/U ratios. In reality, however, the
#' assumption of initial secular equilibrium is rarely met. Accounting
#' for disequilibrium requires a more complex set of age equations,
#' which are based on a coupled system of differetial equations. The
#' solution to this system of equations is given by a matrix
#' exponential. \code{IsoplotR} solves this matrix exponential for any
#' given time, using either the assumed initial activity ratios, or
#' (for young samples) the measured activity ratios of the longest
#' lived intermediate daughters. Based on a \code{Matlab} script by
#' Noah McLean.
#'
#' @param tt the age of the sample
#' @param d an object of class \link{diseq}
#' @param exterr propagate the uncertainties associated with decay
#'     constants and the \eqn{{}^{238}}U/\eqn{{}^{235}}U-ratio.
#'
#' @return a list containing the initial and present-day atomic
#'     abundances of the \eqn{{}^{238}}U-\eqn{{}^{206}}Pb and
#'     \eqn{{}^{235}}U-\eqn{{}^{207}}Pb decay chains; the
#'     \eqn{{}^{206}}Pb/\eqn{{}^{238}}U,
#'     \eqn{{}^{207}}Pb/\eqn{{}^{235}}U and
#'     \eqn{{}^{207}}Pb/\eqn{{}^{206}}Pb ratios at time \code{tt}; the
#'     derivatives of the \eqn{{}^{206}}Pb/\eqn{{}^{238}}U,
#'     \eqn{{}^{207}}Pb/\eqn{{}^{235}}U and
#'     \eqn{{}^{207}}Pb/\eqn{{}^{206}}Pb ratios with respect to time;
#'     and the derivatives of the \eqn{{}^{206}}Pb/\eqn{{}^{238}}U,
#'     \eqn{{}^{207}}Pb/\eqn{{}^{235}}U and
#'     \eqn{{}^{207}}Pb/\eqn{{}^{206}}Pb ratios with respect to the
#'     intermediate decay constants and
#'     \eqn{{}^{238}}U/\eqn{{}^{235}}U-ratio.
#'
#' @seealso \code{\link{diseq}}
#' @author Noah McLean (algorithm) and Pieter Vermeesch (code)
#' @examples
#' d <- diseq(U48=list(x=0,option=1),ThU=list(x=2,option=1),
#'            RaU=list(x=2,option=1),PaU=list(x=2,option=1))
#' mclean(tt=2,d=d)
#' @export
mclean <- function(tt=0,d=diseq(),exterr=FALSE){
    out <- list(truncated=FALSE)
    l38 <- d$L['U238']
    l34 <- d$L['U234']
    l30 <- d$L['Th230']
    l26 <- d$L['Ra226']
    l35 <- d$L['U235']
    l31 <- d$L['Pa231']
    l32 <- lambda('Th232')[1]
    U <- iratio('U238U235')[1]
    nc <- length(d$ThU$x)
    rep(0,nc) -> out$dPb206U238dl38 -> out$dPb206U238dl34 ->
        out$dPb206U238dl30 -> out$dPb206U238dl26 ->
        out$dPb207U235dl35 -> out$dPb207U235dl31
    out$Pb208Th232 <- rep(exp(l32*tt)-1,nc)
    out$dPb208Th232dt <- rep(exp(l32*tt)*l32,nc)
    if (d$equilibrium){
        out$Pb206U238 <- rep(exp(l38*tt)-1,nc)
        out$Pb207U235 <- rep(exp(l35*tt)-1,nc)
        out$dPb206U238dt <- rep(exp(l38*tt)*l38,nc)
        out$dPb207U235dt <- rep(exp(l35*tt)*l35,nc)
        if (exterr){
            out$dPb206U238dl38 <- rep(tt*exp(l38*tt),nc)
            out$dPb207U235dl35 <- rep(tt*exp(l35*tt),nc)
        }
    } else {
        nt <- meas_diseq_nt(d=d,nc=nc)[c('U238','U234','Th230','U235'),,drop=FALSE]
        if (d$U48$option==2){
            t234cutoff <- meas_diseq_maxt(d=d,nuclide='U234')
            n0 <- reverse(tt=min(tt,t234cutoff),mexp=mexp_8405(),nt=nt)
            U48i <- (n0['U234',]*l34)/(n0['U238',]*l38)
            if (any(U48i<0)){
                d$U48 <- list(x=0,sx=0,option=1)
                out$truncated <- TRUE
            }
            if (any(U48i>d$U48$M)){
                d$U48 <- list(x=d$U48$M,sx=0,option=1)
                out$truncated <- TRUE
            }
        }
        if (d$ThU$option==2){
            t230cutoff <- meas_diseq_maxt(d=d,nuclide='Th230')
            n0 <- reverse(tt=min(tt,t230cutoff),mexp=mexp_8405(),nt=nt)
            ThUi <- (n0['Th230',]*l30)/(n0['U238',]*l38)
            if (any(ThUi<0)){
                d$ThU <- list(x=0,sx=0,option=1)
                out$truncated <- TRUE
            }
            if (any(ThUi>d$ThU$M)){
                d$ThU <- list(x=d$ThU$M,sx=0,option=1)
                out$truncated <- TRUE
            }
        }
        n0 <- c(1/l38,1/l34,1/l30,1/l26,0,1/l35,1/l31,0)
        d$n0 <- (n0 %*% matrix(1,nrow=1,ncol=nc)) # duplicate columns
        rownames(d$n0) <- names(d$L)
        if (d$U48$option<2){      # initial 234U
            if (d$U48$option==1) d$n0['U234',] <- d$U48$x/l34
            if (d$ThU$option==1){ # initial 230Th
                d$n0['Th230',] <- d$ThU$x/l30
            } else if (d$ThU$option==2){ # measured 230Th
                nt <- forward(tt=tt,d=d)[c('U238','U234','Th230','U235'),,drop=FALSE]
                nt['Th230',] <- d$ThU$x*nt['U238',]*l38/l30 # overwrite
                d$n0['Th230',] <- reverse(tt=tt,mexp=mexp_8405(),nt=nt)['Th230',]
            }
        } else {                    # measured 234U
            if (d$ThU$option<2){ # initial 230Th
                nt <- forward(tt=tt,d=d)[c('U238','U234','U235'),,drop=FALSE]
                nt['U234',] <- d$U48$x*nt['U238',]*l38/l34 # overwrite
                d$n0['U234',] <- reverse(tt=tt,mexp=mexp_845(),nt=nt)['U234',]
                if (d$ThU$option==1) d$n0['Th230',] <- d$ThU$x/l30
            } else {                # measured 230Th
                nt <- forward(tt=tt,d=d)[c('U238','U234','Th230','U235'),,drop=FALSE]
                nt['U234',] <- d$U48$x*nt['U238',]*l38/l34  # overwrite
                nt['Th230',] <- d$ThU$x*nt['U238',]*l38/l30 # overwrite
                d$n0[c('U234','Th230'),] <-
                    reverse(tt=tt,mexp=mexp_8405(),nt=nt)[c('U234','Th230'),]
            }
        }
        if (d$RaU$option>0) d$n0['Ra226',] <- d$RaU$x/l26
        if (d$PaU$option>0) d$n0['Pa231',] <- d$PaU$x/l31
        d$nt <- forward(tt=tt,d=d)
        dntdt <- forward(tt,d=d,derivative=1)
        out$Pb206U238 <- d$nt['Pb206',]/d$nt['U238',]
        out$Pb207U235 <- d$nt['Pb207',]/d$nt['U235',]
        out$dPb206U238dt <- dntdt['Pb206',]/d$nt['U238',] -
            out$Pb206U238*dntdt['U238',]/d$nt['U238',]
        out$dPb207U235dt <- dntdt['Pb207',]/d$nt['U235',] -
            out$Pb207U235*dntdt['U235',]/d$nt['U235',]
        out$U48i <- (d$n0['U234',]*l34)/(d$n0['U238',]*l38)
        out$ThUi <- (d$n0['Th230',]*l30)/(d$n0['U238',]*l38)
        out$RaUi <- (d$n0['Ra226',]*l26)/(d$n0['U238',]*l38)
        out$PaUi <- (d$n0['Pa231',]*l31)/(d$n0['U235',]*l35)
        out$U48 <- (d$nt['U234',]*l34)/(d$nt['U238',]*l38)
        out$ThU <- (d$nt['Th230',]*l30)/(d$nt['U238',]*l38)
        out$RaU <- (d$nt['Ra226',]*l26)/(d$nt['U238',]*l38)
        out$PaU <- (d$nt['Pa231',]*l31)/(d$nt['U235',]*l35)
        out$n0 <- linkUseries(n=d$n0,U=U*exp((l38-l35)*tt))
        out$nt <- linkUseries(n=d$nt,U=U)
        if (exterr){
            K <- get_diseq_K(tt=tt,d=d)
            out$dPb206U238dl38 <- drdl(d=d,K=K,den='U238',num='Pb206',parent='U238')
            out$dPb206U238dl34 <- drdl(d=d,K=K,den='U238',num='Pb206',parent='U234')
            out$dPb206U238dl30 <- drdl(d=d,K=K,den='U238',num='Pb206',parent='Th230')
            out$dPb206U238dl26 <- drdl(d=d,K=K,den='U238',num='Pb206',parent='Ra226')
            out$dPb207U235dl35 <- drdl(d=d,K=K,den='U235',num='Pb207',parent='U235')
            out$dPb207U235dl31 <- drdl(d=d,K=K,den='U235',num='Pb207',parent='Pa231')
        }
    }
    if (tt>0){
        out$Pb207Pb206 <- out$Pb207U235/(U*out$Pb206U238)
        out$dPb207Pb206dt <- (out$dPb207U235dt*out$Pb206U238 -
                              out$dPb206U238dt*out$Pb207U235)/(U*out$Pb206U238^2)
    } else {
        0 -> out$Pb207Pb206 -> out$dPb207Pb206dt
    }
    if (exterr){
        out$dPb207Pb206dl38 <- -out$dPb206U238dl38*out$Pb207U235/(U*out$Pb206U238^2)
        out$dPb207Pb206dl34 <- -out$dPb206U238dl34*out$Pb207U235/(U*out$Pb206U238^2)
        out$dPb207Pb206dl30 <- -out$dPb206U238dl30*out$Pb207U235/(U*out$Pb206U238^2)
        out$dPb207Pb206dl26 <- -out$dPb206U238dl26*out$Pb207U235/(U*out$Pb206U238^2)
        out$dPb207Pb206dl35 <- out$dPb207U235dl35/(U*out$Pb206U238)
        out$dPb207Pb206dl31 <- out$dPb207U235dl31/(U*out$Pb206U238)
        out$dPb208Th232dl32 <- exp(l32*tt)*tt
        out$dPb207Pb206dU <- -out$Pb207Pb206/U
    } else {
        rep(0,nc) -> out$dPb207Pb206dl38 -> out$dPb207Pb206dl34 ->
            out$dPb207Pb206dl30 -> out$dPb207Pb206dl26 ->
            out$dPb207Pb206dl35 -> out$dPb207Pb206dl31 ->
            out$dPb208Th232dl32 -> out$dPb207Pb206dU
    }
    out
}

linkUseries <- function(n,U){
    num <- n[c('U238','U234','Th230','Ra226','Pb206'),,drop=FALSE]
    den <- n['U238',,drop=FALSE]
    U8series <- sweep(num,2,den,'/')
    num <- n[c('U235','Pa231','Pb207'),,drop=FALSE]
    den <- n['U235',,drop=FALSE]
    U5series <- sweep(num,2,den,'/')
    rbind(U*U8series,U5series)
}

get_diseq_K <- function(tt=0,d=diseq()){
    nl <- length(d$L)
    out <- matrix(0,nl,nl)
    colnames(out) <- names(d$L)
    rownames(out) <- names(d$L)
    for (i in 1:nl){
        for (j in 1:nl){
            if (i==j){
                out[i,j] <- tt*exp(-d$L[i]*tt)
            } else if (d$L[i]==d$L[j]){
                # do nothing
            } else {
                out[i,j] <- (exp(-d$L[i]*tt)-exp(-d$L[j]*tt))/(d$L[j]-d$L[i])
            }
        }
    }
    out
}
drdl <- function(tt=0,K=matrix(0,8,8),d=diseq(),
                 num='Pb206',den='U238',parent='Th230'){
    # 1. create matrix derivative
    i <- which(names(d$L)%in%parent)
    nl <- length(d$L)
    dAdl <- matrix(0,nl,nl)
    dAdl[c(i,i+1),i] <- c(-1,1)
    H <- d$Qinv %*% dAdl %*% d$Q
    P <- H * K
    dntdl <- (d$Q %*% P %*% d$Qinv %*% d$nt)
    rownames(dntdl) <- names(d$L)
    # 2. derivative of the ratio
    out <- (d$nt[den,]*dntdl[num,]-d$nt[num,]*dntdl[den,])/d$nt[den,,drop=FALSE]^2
    as.vector(out)
}
diseq_75_misfit <- function(tt,x,d){
    pred <- subset(age_to_Pb207U235_ratio(tt,d=d),select='75')
    get.5678.misfit(obs=x,pred)
}
diseq_68_misfit <- function(tt,x,d){
    pred <- subset(age_to_Pb206U238_ratio(tt,d=d),select='68')
    get.5678.misfit(obs=x,pred)
}
get.76.misfit <- function(tt,x,d=diseq()){
    pred <- subset(age_to_Pb207Pb206_ratio(tt=tt,d=d),select='76')
    get.5678.misfit(obs=x,pred)
}
get.5678.misfit <- function(obs,pred){
    (log(obs)-log(pred))^2
}

meas_diseq_maxt <- function(d,nuclide='auto'){
    if ((d$ThU$option==2 && nuclide!='U234') || nuclide=='Th230'){
        out <- 20*log(2)/d$L['Th230']
    } else if (d$U48$option==2 || nuclide=='U234'){
        out <- 20*log(2)/d$L['U234']
    } else {
        out <- 4500
    }
    out
}
concordia_end <- function(d,nuclide='auto'){
    tlim <- c(0,meas_diseq_maxt(d=d,nuclide=nuclide))
    nt <- meas_diseq_nt(d=d)[c('U238','U234','Th230','U235'),,drop=FALSE]
    if ((d$ThU$option==2 && nuclide!='U234') || nuclide=='Th230'){
        misfit <- function(tt,d){
            n0 <- reverse(tt,mexp=mexp_8405(),nt=nt)
            n0['Th230',]^2
        }
        out <- stats::optimise(misfit,tlim,d=d)$minimum
    } else if (d$U48$option==2 || nuclide=='U234'){
        misfit <- function(tt,d){
            n0 <- reverse(tt,mexp=mexp_8405(),nt=nt)
            U48i <- (n0['U234',]*d$L['U234'])/(n0['U238',]*d$L['U238'])
            ifelse(d$U48$x<1,U48i^2,(U48i-20)^2)
        }
        out <- stats::optimise(misfit,tlim,d=d)$minimum
    } else {
        out <- 4500
    }
    out
}

measured2initial <- function(x,fit){
    out <- x
    fnames <- names(fit$par)
    dnames <- c('U48','ThU','RaU','PaU')
    for (dname in dnames){
        iname <- paste0(dname,'i')
        if (iname %in% fnames){
            out$d[[dname]]$x <- fit$par[iname]
            out$d[[dname]]$sx <- sqrt(fit$cov[iname,iname])
            out$d[[dname]]$option <- 1
        }
    }
    out
}

meas_diseq_nt <- function(d,nc=1){
    nt <- 0*d$L
    nt[d$L>0] <- 1/d$L[d$L>0]
    if (d$U48$option==2){
        nt['U234'] <- d$U48$x*nt['U238']*d$L['U238']/d$L['U234']
    }
    if (d$ThU$option==2){
        nt['Th230'] <- d$ThU$x*nt['U238']*d$L['U238']/d$L['Th230']
    }
    out <- (nt %*% matrix(1,nrow=1,ncol=nc)) # duplicate columns
    rownames(out) <- names(d$L)
    out
}

pname2aname <- function(pname){
    substr(pname, 1, nchar(pname)-1)
}
aname2pname <- function(aname){
    paste0(aname,'i')
}
