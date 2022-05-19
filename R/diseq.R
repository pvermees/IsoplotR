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
#' @param U48 a list containing three items (\code{x}, \code{sx} and
#'     \code{option}) specifying the \eqn{{}^{234}}U/\eqn{{}^{238}}U
#'     disequilibrium.
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
#' @param ThU a list containing three items (\code{x}, \code{sx} and
#'     \code{option}) specifying the \eqn{{}^{230}}Th/\eqn{{}^{238}}U
#'     disequilibrium.
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
#' @param RaU a list containing three items (\code{x}, \code{sx} and
#'     \code{option}) specifying the \eqn{{}^{226}}Ra/\eqn{{}^{238}}U
#'     disequilibrium.
#'
#' If \code{option=0}, then \code{x} and \code{sx} are ignored and no
#'     disequilibrium correction is applied.
#'
#' If \code{option=1}, then \code{x} contains the initial
#'     \eqn{{}^{226}}Ra/\eqn{{}^{238}}U ratio and \code{sx} its
#'     standard error.
#' 
#' @param PaU a list containing three items (\code{x}, \code{sx} and
#'     \code{option}) specifying the \eqn{{}^{231}}Pa/\eqn{{}^{235}}U
#'     disequilibrium.
#'
#' If \code{option=0}, then \code{x} and \code{sx} are ignored and no
#'     disequilibrium correction is applied.
#'
#' If \code{option=1}, then \code{x} contains the initial
#'     \eqn{{}^{231}}Pa/\eqn{{}^{235}}U ratio and \code{sx} its
#'     standard error.
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
#' \item{n0}{the initial atomic abundances of all the parent and
#' daughter isotopes (used by \code{\link{mclean}})}
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
#' concordia(UPb,type=2,xlim=c(0,5000),ylim=c(0.047,0.057))
#' @export
diseq <- function(U48=list(x=1,sx=0,option=0),
                  ThU=list(x=1,sx=0,option=0),
                  RaU=list(x=1,sx=0,option=0),
                  PaU=list(x=1,sx=0,option=0)){
    out <- list()
    class(out) <- 'diseq'
    out$U48 <- U48
    out$ThU <- ThU
    out$RaU <- RaU
    out$PaU <- PaU
    out$equilibrium <- check.equilibrium(d=out)
    l38 <- settings('lambda','U238')[1]
    l34 <- settings('lambda','U234')[1]*1000
    l30 <- settings('lambda','Th230')[1]*1000
    l26 <- settings('lambda','Ra226')[1]*1000
    l35 <- settings('lambda','U235')[1]
    l31 <- settings('lambda','Pa231')[1]*1000
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
    out$L <- c(l38,l34,l30,l26,0,l35,l31,0)
    out$n0 <- matrix(c(1/l38,1/l34,1/l30,1/l26,0,1/l35,1/l31,0),ncol=1)
    nuclides <- c('U238','U234','Th230','Ra226','Pb206','U235','Pa231','Pb207')
    names(out$L) <- nuclides
    rownames(out$n0) <- nuclides
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

geomean.diseq <- function(x,...){
    out <- x
    for (ratio in c('U48','ThU','RaU','PaU')){
        out[[ratio]]$x <- geomean(x[[ratio]]$x,...)
    }
    out
}

# to infer initial 38,34,30,35 from measured 34/38 and 30/38
mexp.8405 <- function(){
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
mexp.845 <- function(nratios=3){
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
    out <- (mexp$Q %*% diag(exp(mexp$L*tt)) %*% mexp$Qinv %*% nt)
    rownames(out) <- rownames(nt)
    out
}
forward <- function(tt,d=diseq(),derivative=0){
    if (derivative==0){
        out <- (d$Q %*% diag(exp(-d$L*tt)) %*% d$Qinv %*% d$n0)
    } else if (derivative==1){
        out <- (d$Q %*% diag(exp(-d$L)) %*%
                diag(exp(-d$L*tt)) %*% d$Qinv %*% d$n0)
    } else if (derivative==2){
        out <- (d$Q %*% diag(exp(-d$L)) %*% diag(exp(-d$L)) %*%
                diag(exp(-d$L*tt)) %*% d$Qinv %*% d$n0)
    }
    rownames(out) <- names(d$L)
    out
}

check.equilibrium <- function(d=diseq()){
    U48 <- (d$U48$option==0 | all(d$U48$x==1))
    ThU <- (d$ThU$option==0 | all(d$ThU$x==1))
    RaU <- (d$RaU$option==0 | all(d$RaU$x==1))
    PaU <- (d$PaU$option==0 | all(d$PaU$x==1))
    U48 & ThU & RaU & PaU
}
measured.disequilibrium <- function(d=diseq()){
    equilibrium <- check.equilibrium(d)
    measured <- (d$U48$option==2 | d$ThU$option==2)
    !equilibrium & measured
}

#' @title Predict disequilibrium concordia compositions
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
#' @return
#' a list containing the predicted \eqn{{}^{206}}Pb/\eqn{{}^{238}}U,
#' \eqn{{}^{207}}Pb/\eqn{{}^{235}}U and
#' \eqn{{}^{207}}Pb/\eqn{{}^{206}}Pb ratios at time \code{tt}; the
#' derivatives of the \eqn{{}^{206}}Pb/\eqn{{}^{238}}U,
#' \eqn{{}^{207}}Pb/\eqn{{}^{235}}U and
#' \eqn{{}^{207}}Pb/\eqn{{}^{206}}Pb ratios with respect to time; and
#' the derivatives of the \eqn{{}^{206}}Pb/\eqn{{}^{238}}U,
#' \eqn{{}^{207}}Pb/\eqn{{}^{235}}U and
#' \eqn{{}^{207}}Pb/\eqn{{}^{206}}Pb ratios with respect to the
#' intermediate decay constants and
#' \eqn{{}^{238}}U/\eqn{{}^{235}}U-ratio.
#'
#' @seealso \code{\link{diseq}}
#' @author Noah McLean (algorithm) and Pieter Vermeesch (code)
#' @examples
#' d <- diseq(U48=list(x=0,option=1),ThU=list(x=2,option=1),
#'            RaU=list(x=2,option=1),PaU=list(x=2,option=1))
#' mclean(tt=2,d=d)
#' @export
mclean <- function(tt=0,d=diseq(),exterr=FALSE){
    out <- list()
    l38 <- lambda('U238')[1]
    l34 <- lambda('U234')[1]*1000
    l30 <- lambda('Th230')[1]*1000
    l26 <- lambda('Ra226')[1]*1000
    l35 <- lambda('U235')[1]
    l32 <- lambda('Th232')[1]
    l31 <- lambda('Pa231')[1]*1000
    U <- iratio('U238U235')[1]
    nc <- length(d$ThU$x)
    out$dPb206U238dl38 <- rep(0,nc)
    out$dPb206U238dl34 <- rep(0,nc)
    out$dPb206U238dl30 <- rep(0,nc)
    out$dPb206U238dl26 <- rep(0,nc)
    out$dPb207U235dl35 <- rep(0,nc)
    out$dPb207U235dl31 <- rep(0,nc)
    out$Pb208Th232 <- rep(exp(l32*tt)-1,nc)
    out$dPb208Th232dt <- rep(exp(l32*tt)*l32,nc)
    out$d2Pb208Th232dt2 <- rep(exp(l32*tt)*l32^2,nc)
    if (check.equilibrium(d=d)){
        out$Pb206U238 <- rep(exp(l38*tt)-1,nc)
        out$Pb207U235 <- rep(exp(l35*tt)-1,nc)
        out$dPb206U238dt <- rep(exp(l38*tt)*l38,nc)
        out$dPb207U235dt <- rep(exp(l35*tt)*l35,nc)
        out$d2Pb206U238dt2 <- rep(exp(l38*tt)*l38^2,nc)
        out$d2Pb207U235dt2 <- rep(exp(l35*tt)*l35^2,nc)
        if (exterr){
            out$dPb206U238dl38 <- rep(tt*exp(l38*tt),nc)
            out$dPb207U235dl35 <- rep(tt*exp(l35*tt),nc)
        }
    } else {
        d$n0 <- (d$n0 %*% matrix(1,nrow=1,ncol=nc)) # duplicate columns
        rownames(d$n0) <- names(d$L)
        if (d$U48$option<2){      # initial 234U
            if (d$U48$option==1) d$n0['U234',] <- d$U48$x/l34
            if (d$ThU$option==1){ # initial 230Th
                d$n0['Th230',] <- d$ThU$x/l30
            } else if (d$ThU$option==2){ # measured 230Th
                nt <- forward(tt=tt,d=d)[c('U238','U234','Th230','U235'),,drop=FALSE]
                nt['Th230',] <- d$ThU$x*nt['U238',]*l38/l30 # overwrite
                d$n0['Th230',] <-
                    reverse(tt=tt,mexp=mexp.8405(),nt=nt)['Th230',]
            }
        } else {                 # measured 234U
            if (d$ThU$option<2){ # initial 230Th
                nt <- forward(tt=tt,d=d)[c('U238','U234','U235'),,drop=FALSE]
                nt['U234',] <- d$U48$x*nt['U238',]*l38/l34 # overwrite
                d$n0['U234',] <- reverse(tt=tt,mexp=mexp.845(),nt=nt)['U234',]
                if (d$ThU$option==1) d$n0['Th230',] <- d$ThU$x/l30
            } else {             # measured 230Th
                nt <- forward(tt=tt,d=d)[c('U238','U234','Th230','U235'),,drop=FALSE]
                nt['U234',] <- d$U48$x*nt['U238',]*l38/l34 # overwrite
                nt['Th230',] <- d$ThU$x*nt['U238',]*l38/l30 # overwrite
                d$n0[c('U234','Th230'),] <-
                    reverse(tt=tt,mexp=mexp.8405(),nt=nt)[c('U234','Th230'),]
            }
        }
        if (d$RaU$option>0) d$n0['Ra226',] <- d$RaU$x/l26
        if (d$PaU$option>0) d$n0['Pa231',] <- d$PaU$x/l31
        out$U48i <- (d$n0['U234',]*l34)/(d$n0['U238',]*l38)
        out$ThUi <- (d$n0['Th230',]*l30)/(d$n0['U238',]*l38)
        d$nt <- forward(tt=tt,d=d)
        dntdt <- forward(tt,d=d,derivative=1)
        d2ntdt2 <- forward(tt,d=d,derivative=2)
        out$Pb206U238 <- d$nt['Pb206',]/d$nt['U238',]
        out$Pb207U235 <- d$nt['Pb207',]/d$nt['U235',]
        out$dPb206U238dt <- dntdt['Pb206',]/d$nt['U238',] -
            out$Pb206U238*dntdt['U238',]/d$nt['U238',]
        out$d2Pb206U238dt2 <- d2ntdt2['Pb206',]/d$nt['U238',] -
            2*dntdt['Pb206',]*dntdt['U238',]/d$nt['U238',]^2 -
            out$Pb206U238*d2ntdt2['U238',]/d$nt['U238',] +
            2*out$Pb206U238*(dntdt['U238',]/d$nt['U238',])^2
        out$dPb207U235dt <- dntdt['Pb207',]/d$nt['U235',] -
            out$Pb207U235*dntdt['U235',]/d$nt['U235',]
        out$d2Pb207U235dt2 <- d2ntdt2['Pb207',]/d$nt['U235',] -
            2*dntdt['Pb207',]*dntdt['U235',]/d$nt['U235',]^2 -
            out$Pb207U235*d2ntdt2['U235',]/d$nt['U235',] +
            2*out$Pb207U235*(dntdt['U235',]/d$nt['U235',])^2
        if (exterr){
            K <- get.diseq.K(tt=tt,d=d)
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
        out$Pb207Pb206 <- rep(0,nc)
        out$dPb207Pb206dt <- rep(0,nc)
    }
    if (exterr){
        out$dPb207Pb206dl38 <- -out$dPb206U238dl38*out$Pb207U235/(U*out$Pb206U238^2)
        out$dPb207Pb206dl34 <- -out$dPb206U238dl34*out$Pb207U235/(U*out$Pb206U238^2)
        out$dPb207Pb206dl30 <- -out$dPb206U238dl30*out$Pb207U235/(U*out$Pb206U238^2)
        out$dPb207Pb206dl26 <- -out$dPb206U238dl26*out$Pb207U235/(U*out$Pb206U238^2)
        out$dPb207Pb206dl35 <- out$dPb207U235dl35*out$Pb206U238/(U*out$Pb206U238^2)
        out$dPb207Pb206dl31 <- out$dPb207U235dl31*out$Pb206U238/(U*out$Pb206U238^2)
        out$dPb208Th232dl32 <- rep(exp(l32*tt)*tt,nc)
        out$dPb207Pb206dU <- -out$Pb207Pb206/U
    } else {
        out$dPb207Pb206dl38 <- rep(0,nc)
        out$dPb207Pb206dl34 <- rep(0,nc)
        out$dPb207Pb206dl30 <- rep(0,nc)
        out$dPb207Pb206dl26 <- rep(0,nc)
        out$dPb207Pb206dl35 <- rep(0,nc)
        out$dPb207Pb206dl31 <- rep(0,nc)
        out$dPb208Th232dl32 <- rep(0,nc)
        out$dPb207Pb206dU <- rep(0,nc)
    }
    out
}

get.diseq.K <- function(tt=0,d=diseq()){
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
    out <- (d$nt[den,]*dntdl[num,]-
            d$nt[num,]*dntdl[den,])/d$nt[den,,drop=FALSE]^2
    out
}

diseq.75.misfit <- function(tt,x,d){
    pred <- subset(age_to_Pb207U235_ratio(tt,d=d),select='75')
    get.5678.misfit(obs=x,pred)
}
diseq.68.misfit <- function(tt,x,d){
    pred <- subset(age_to_Pb206U238_ratio(tt,d=d),select='68')
    get.5678.misfit(obs=x,pred)
}
get.76.misfit <- function(tt,x,d=diseq()){
    pred <- subset(age_to_Pb207Pb206_ratio(tt=tt,d=d),select='76')
    get.5678.misfit(obs=x,pred)
}
get.5678.misfit <- function(obs,pred){
    abs(log(obs+1) - log(pred+1)) # +1 for robustness to negative values
}

meas.diseq.maxt <- function(d){
    if (d$ThU$option==2){
        out <- stats::optimise(function(tt,d) mclean(tt=tt,d=d)$ThUi^2,
                               c(0,1),d=d)$minimum
    } else if (d$U48$option==2){
        M <- ifelse(d$U48$x<1,0,500)
        out <- stats::optimise(function(tt,d,maxU48) (mclean(tt=tt,d=d)$U48i-maxU48)^2,
                               c(0,10),d=d,maxU48=M)$minimum
    } else {
        out <- 4500
    }
    out
}
