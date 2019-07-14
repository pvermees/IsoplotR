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
#' between the activity of \eqn{{}^{238}}U, \eqn{{}^{234}}Th,
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
#'     \eqn{{}^{234}}U/\eqn{{}^{238}}U ratio.
#'
#' If \code{option=2}, then \code{x} contains the measured
#'     \eqn{{}^{234}}U/\eqn{{}^{238}}U ratio.
#'
#' @param ThU a list containing three items (\code{x}, \code{sx} and
#'     \code{option}) specifying the \eqn{{}^{230}}Th/\eqn{{}^{238}}U
#'     disequilibrium.
#'
#' If \code{option=0}, then \code{x} and \code{sx} are ignored and no
#'     disequilibrium correction is applied.
#'
#' If \code{option=1}, then \code{x} contains the initial
#'     \eqn{{}^{230}}Th/\eqn{{}^{238}}U ratio.
#'
#' If \code{option=2}, then \code{x} contains the measured
#'     \eqn{{}^{230}}Th/\eqn{{}^{238}}U ratio.
#'
#' If \code{option=3}, then \code{x} contains the measured Th/U ratio
#'     of the magma (assumed or determined from the whole rock or
#'     volcanic glass). This only applies for Th-bearing U-Pb data
#'     formats 7 and 8.
#'
#' @param RaU a list containing three items (\code{x}, \code{sx} and
#'     \code{option}) specifying the \eqn{{}^{226}}Ra/\eqn{{}^{238}}U
#'     disequilibrium.
#'
#' If \code{option=0}, then \code{x} and \code{sx} are ignored and no
#'     disequilibrium correction is applied.
#'
#' If \code{option=1}, then \code{x} contains the initial
#'     \eqn{{}^{226}}Ra/\eqn{{}^{238}}U ratio.
#'
#' If \code{option=2}, then \code{x} contains the measured
#'     \eqn{{}^{226}}Ra/\eqn{{}^{238}}U ratio.
#' 
#' @param PaU a list containing three items (\code{x}, \code{sx} and
#'     \code{option}) specifying the \eqn{{}^{231}}Pa/\eqn{{}^{235}}U
#'     disequilibrium.
#'
#' If \code{option=0}, then \code{x} and \code{sx} are ignored and no
#'     disequilibrium correction is applied.
#'
#' If \code{option=1}, then \code{x} contains the initial
#'     \eqn{{}^{231}}Pa/\eqn{{}^{235}}U ratio.
#'
#' If \code{option=2}, then \code{x} contains the measured
#'     \eqn{{}^{231}}Ra/\eqn{{}^{235}}U ratio.
#' 
#' @return a list with the activity ratios, an eigen composition of
#'     the decay contant matrix and the atomic abundances of the
#'     parent and (intermediate) daughter nuclides
#' 
#' @examples
#' d <- diseq(U48=list(x=0,option=1),ThU=list(x=2,option=1),
#'            RaU=list(x=2,option=1),PaU=list(x=2,option=1))
#' fn <- system.file("diseq.csv",package="IsoplotR")
#' UPb <- read.data(fn,method='U-Pb',format=2)
#' concordia(UPb,type=2,show.age=1)
#' @export
diseq <- function(U48=list(x=1,sx=0,option=0),
                  ThU=list(x=1,sx=0,option=0),
                  RaU=list(x=1,sx=0,option=0),
                  PaU=list(x=1,sx=0,option=0)){
    out <- list()
    class(out) <- 'diseq'
    out$equilibrium <- (U48$option==0 & ThU$option==0 & RaU$option==0 & PaU$option==0) |
        (U48$x==1 & ThU$x==0 & RaU$x==0 & PaU$x==0)
    out$U48 <- U48
    out$ThU <- ThU
    out$RaU <- RaU
    out$PaU <- PaU
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
    out$Qinv[8,6:8] <-1
    nuclides <- c('U238','U234','Th230','Ra226','Pb206','U235','Pa231','Pb207')
    out$L <- c(l38,l34,l30,l26,0,l35,l31,0)
    out$n0 <- rep(0,8)
    out$nt <- rep(0,8)
    names(out$n0) <- nuclides
    names(out$nt) <- nuclides
    names(out$L) <- nuclides
    out
}

#' @export
`[.diseq` <- function(x,i){
    out <- x
    for (ratio in c('U48','ThU','RaU','PaU')){
        j <- min(length(out[[ratio]]$x),i)
        out[[ratio]]$x <- out[[ratio]]$x[j]
    }
    out
}

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
        out[[ratio]]$x <- geomean(x[[ratio]]$x)
    }
    out
}

get_atoms <- function(tt=0,d=diseq(),nuclide='Th230',ratio='ThU',parent=FALSE){
    out <- d
    if (parent) init <- 1
    else init <- d[[ratio]]$x
    if (d[[ratio]]$option==0){
        out$n0[nuclide] <- 1/d$L[nuclide]
        out$nt[nuclide] <- forward(tt=tt,d=out)[nuclide]
    } else if (d[[ratio]]$option==1){
        out$n0[nuclide] <- init/d$L[nuclide]
        out$nt[nuclide] <- forward(tt=tt,d=out)[nuclide]
    } else if (tt*d$L[nuclide]<10){
        out$nt[nuclide] <- init/d$L[nuclide]
        out$n0[nuclide] <- reverse(tt=tt,d=out)[nuclide]        
    } else { # measured and expired
        out$n0[nuclide] <- 1/d$L[nuclide]
        out$nt[nuclide] <- forward(tt=tt,d=out)[nuclide]        
    }
    out
}

# d only contains one sample
mclean <- function(tt=0,d=diseq()){
    out <- list()
    if (d$equilibrium){
        out$Pb206U238 <- exp(d$L['U238']*tt) - 1
        out$Pb207U235 <- exp(d$L['U235']*tt) - 1
        out$dPb206U238dt <- exp(d$L['U238']*tt)*d$L['U238']
        out$dPb207U235dt <- exp(d$L['U235']*tt)*d$L['U235']
        out$d2Pb206U238dt2 <- exp(d$L['U238']*tt)*d$L['U238']^2
        out$d2Pb207U235dt2 <- exp(d$L['U235']*tt)*d$L['U235']^2
        out$dPb206U238dl8 <- tt*exp(d$L['U238']*tt)
        out$dPb207U235dl5 <- tt*exp(d$L['U235']*tt)
    } else {
        d <- get_atoms(tt=tt,d=d,nuclide='U238',ratio='U48',parent=TRUE)
        d <- get_atoms(tt=tt,d=d,nuclide='U234',ratio='U48')
        d <- get_atoms(tt=tt,d=d,nuclide='Th230',ratio='ThU')
        d <- get_atoms(tt=tt,d=d,nuclide='Ra226',ratio='RaU')
        d <- get_atoms(tt=tt,d=d,nuclide='U235',ratio='PaU',parent=TRUE)
        d <- get_atoms(tt=tt,d=d,nuclide='Pa231',ratio='PaU')
        nt <- forward(tt=tt,d=d)
        dntdt <- forward(tt,d=d,derivative=1)
        d2ntdt2 <- forward(tt,d=d,derivative=2)
        out$Pb206U238 <- nt['Pb206']/nt['U238']
        out$Pb207U235 <- nt['Pb207']/nt['U235']
        out$dPb206U238dt <- dntdt['Pb206']/nt['U238'] -
            out$Pb206U238*dntdt['U238']/nt['U238']
        out$d2Pb206U238dt2 <- d2ntdt2['Pb206']/nt['U238'] -
            2*dntdt['Pb206']*dntdt['U238']/nt['U238']^2 -
            out$Pb206U238*d2ntdt2['U238']/nt['U238'] +
            2*out$Pb206U238*(dntdt['U238']/nt['U238'])^2
        out$dPb207U235dt <- dntdt['Pb207']/nt['U235'] -
            out$Pb207U235*dntdt['U235']/nt['U235']
        out$d2Pb207U235dt2 <- d2ntdt2['Pb207']/nt['U235'] -
            2*dntdt['Pb207']*dntdt['U235']/nt['U235']^2 -
            out$Pb207U235*d2ntdt2['U235']/nt['U235'] +
            2*out$Pb207U235*(dntdt['U235']/nt['U235'])^2
        out$dPb206U238dl8 <- 0 # TODO
        out$dPb207U235dl5 <- 0 # TODO
    }
    out
}
forward <- function(tt,d=diseq(),derivative=0){
    if (derivative==0){
        out <- as.vector(d$Q %*% diag(exp(-d$L*tt)) %*% d$Qinv %*% d$n0)
    } else if (derivative==1){
        out <- as.vector(d$Q %*% diag(exp(-d$L)) %*%
                         diag(exp(-d$L*tt)) %*% d$Qinv %*% d$n0)
    } else if (derivative==2){
        out <- as.vector(d$Q %*% diag(exp(-d$L)) %*% diag(exp(-d$L)) %*%
                         diag(exp(-d$L*tt)) %*% d$Qinv %*% d$n0)
    }
    names(out) <- names(d$n0)
    out
}
reverse <- function(tt,d=diseq(),derivative=FALSE){
    if (derivative==0){
        out <- as.vector(d$Q %*% diag(exp(d$L*tt)) %*% d$Qinv %*% d$nt)
    } else if (derivative==1){
        out <- as.vector(d$Q %*% diag(exp(d$L)) %*%
                         diag(exp(d$L*tt)) %*% d$Qinv %*% d$nt)
    } else if (derivative==2){
        out <- as.vector(d$Q %*% diag(exp(d$L)) %*% diag(exp(d$L)) %*%
                         diag(exp(d$L*tt)) %*% d$Qinv %*% d$nt)        
    }
    names(out) <- names(d$nt)
    out
}

diseq.75.misfit <- function(tt,x,d){
    abs(log(x) - log(subset(age_to_Pb207U235_ratio(tt,d=d),select='75')))
}
diseq.68.misfit <- function(tt,x,d){
    abs(log(x) - log(subset(age_to_Pb206U238_ratio(tt,d=d),select='68')))
}
Pb207Pb206.misfit <- function(tt,x,d=diseq()){
    abs(x - subset(age_to_Pb207Pb206_ratio(tt,d=d),select='76'))
}
# derivative of the 7/6 misfit function w.r.t. time
dmf76dt <- function(x,t.76,d=diseq()){
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    U <- iratio('U238U235')[1]
    D <- mclean(tt=t.76,d=d)
    r75 <- D$Pb207U235
    r68 <- D$Pb206U238
    dr75dt <- D$dPb207U235dt
    dr68dt <- D$dPb206U238dt
    2*(r75/(U*r68)-x)*(dr75dt*r68-r75*dr68dt)/(U*r68^2)
}
dmf76dl5 <- function(x,t.76,d=diseq()){
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    U <- iratio('U238U235')[1]
    D <- mclean(tt=t.76,d=d)
    r75 <- D$Pb207U235
    r68 <- D$Pb206U238
    dr75dl5 <- D$dPb207U235dl5
    2*(r75/(U*r68)-x)/(U*r68)
}
dmf76dl8 <- function(x,t.76,d=diseq()){
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    U <- iratio('U238U235')[1]
    D <- mclean(tt=t.76,d=d)
    r75 <- D$Pb207U235
    r68 <- D$Pb206U238
    dr68dl8 <- D$dPb206U238dl8
    2*(r75/(U*r68)-x)*(-r75*dr68dl8)/(U*r68^2)
}
dmf76dU <- function(x,t.76,d=diseq()){
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    U <- iratio('U238U235')[1]
    D <- mclean(tt=t.76,d=d)
    r75 <- D$Pb207U235
    r68 <- D$Pb206U238
    -2*(r75/(U*r68)-x)*r75/(r68*U^2)
}
