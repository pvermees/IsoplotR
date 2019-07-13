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
#' There are four ways to correct for the initial disequilibrium
#' between the activity of \eqn{{}^{238}}U, \eqn{{}^{234}}Th,
#' \eqn{{}^{230}}Th, and \eqn{{}^{226}}Ra; or between \eqn{{}^{235}}U
#' and \eqn{{}^{231}}Pa:
#'
#' \enumerate{
#' 
#' \item{Specify the assumed initial activity ratios and calculate how
#' much excess \eqn{{}^{206}}Pb and \eqn{{}^{207}}Pb these would have
#' produced. (Wendt and Carl, 1985).}
#' 
#' \item{Measure the current activity ratios to infer the initial
#' ratios.  This approach only works for young samples (< 5Ma, say).}
#' 
#' \item{Specify the elemental fractionation factor between Th and U
#' in the magma chamber and the mineral (Schaerer,
#' 1984). \code{IsoplotR} generalises this approach to Ra/U and Pa/U
#' as well. However, it still assumes secular equilibrium between
#' \eqn{{}^{234}}U and \eqn{{}^{238}}Th.}
#' 
#' }
#' 
#' @param option one of four options:
#'
#' \code{0}: no disequilibrium correction
#'
#' \code{1}: use assumed initial activity ratios
#'
#' \code{2}: use measured current activity ratios
#'
#' \code{3}: use partition coefficients between the mineral and magma
#'
#' \code{4}: use the measured Th/U ratio of the magma to compute the
#' Th/U partition coefficient. This option only works for U-Pb
#' datasets of format 7 or 8.
#'
#' @param U48 the \eqn{^{234}}U/\eqn{^{238}}U-activity ratio (initial
#'     if \code{option=1} or measured if \code{option=2}).
#' 
#' @param Th0U8 the \eqn{^{230}}Th/\eqn{^{238}}U-activity ratio
#'     (initial if \code{option=1} or measured if \code{option=2}).
#' 
#' @param Ra6U8 the \eqn{^{226}}Ra/\eqn{^{238}}U-activity ratio
#'     (initial if \code{option=1} or measured if \code{option=2}).
#' 
#' @param Pa1U5 the \eqn{^{231}}Pa/\eqn{^{235}}U-activity ratio
#'     (initial if \code{option=1} or measured if \code{option=2}).
#' 
#' @param fThU the Th/U fractionation factor between the mineral (m)
#'     and the magma (M): \code{fThU} = (Th/U)\eqn{_m}/(Th/U)\eqn{_M}.
#' 
#' @param fRaU the Ra/U fractionation factor between the mineral (m)
#'     and the magma (M): \code{fRaU} = (Ra/U)\eqn{_m}/(Ra/U)\eqn{_M}.
#' 
#' @param fPaU the Pa/U fractionation factor between the mineral (m)
#'     and the magma (M): \code{fPaU} = (Pa/U)\eqn{_m}/(Pa/U)\eqn{_M}.
#'  
#' @param ThU the measured Th/U ratio of the magma, which is to be
#'     combined with the Th/U ratios of the minerals to compute the
#'     Th/U fractionation factor.
#'  
#' @return a list with the following items: \code{option} and
#'     (\code{U48}, \code{Th08}, \code{Ra6U8}, \code{Pa1U8}) [if
#'     \code{option=1} or \code{option=2}] and (\code{fThU},
#'     \code{RaU}, \code{PaU}) [if \code{option=3}].
#' 
#' @examples
#' d <- diseq(option=3,fThU=2)
#' fn <- system.file("UPb1.csv",package="IsoplotR")
#' UPb <- read.data(fn,method='U-Pb',format=1,d=d)
#' concordia(UPb)
#' 
#' @references
#' Schaerer, U., 1984. The effect of initial \eqn{{}^{230}}Th
#' disequilibrium on young UPb ages: the Makalu case, Himalaya. Earth
#' and Planetary Science Letters, 67(2), pp.191-204.
#' 
#' Wendt, I. and Carl, C., 1985. U/Pb dating of discordant 0.1 Ma old
#' secondary U minerals. Earth and Planetary Science Letters, 73(2-4),
#' pp.278-284.
#' @export
diseq <- function(U48=list(x=1,sx=0,option=0),
                  ThU=list(x=1,sx=0,option=0),
                  RaU=list(x=1,sx=0,option=0),
                  PaU=list(x=1,sx=0,option=0)){
    out <- list()
    class(out) <- 'diseq'
    out$equilibrium <- U48$option==0 & ThU$option==0 & RaU$option==0 & PaU$option==0
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
    out$L <- c(l38,l34,l30,l26,0,l35,l31,0)
    out$atoms <- rep(0,8) # initialise
    names(out$atoms) <- c('U238','U234','Th230','Ra226',
                          'Pb206','U235','Pa231','Pb207')
    names(out$L) <- names(out$atoms)
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
        d$atoms[1] <- 1/d$L['U238']
        d$atoms[2] <- d$U48[[1]]/d$L['U234']
        d$atoms[3] <- d$ThU[[1]]/d$L['Th230']
        d$atoms[4] <- d$RaU[[1]]/d$L['Ra226']
        d$atoms[6] <- 1/d$L['U235']
        d$atoms[7] <- d$PaU[[1]]/d$L['Pa231']
        factor <- 20
        if (d$U48$option>1 & tt*d$L['U234']<factor)
            d$atoms['U234'] <- reverse(tt=tt,d=d)['U234']
        if (d$ThU$option>1 & tt*d$L['Th230']<factor)
            d$atoms['Th230'] <- reverse(tt=tt,d=d)['Th230']
        if (d$RaU$option>1 & tt*d$L['Ra226']<factor)
            d$atoms['Ra226'] <- reverse(tt=tt,d=d)['Ra226']
        if (d$PaU$option>1 & tt*d$L['Pa231']<factor)
            d$atoms['Pa231'] <- reverse(tt=tt,d=d)['Pa231']
        nt <- forward(tt,d=d)
        dntdt <- forward(tt,d=d,derivative=1)
        d2ntdt2 <- forward(tt,d=d,derivative=2)
        out$Pb206U238 <- nt['Pb206']/nt['U238']
        out$Pb207U235 <- nt['Pb207']/nt['U235']
        out$dPb206U238dt <- dntdt['Pb206']/nt['U238'] -
            out$Pb206U238*dntdt['U238']/nt['U238']
        out$d2Pb206U238dt2 <- d2ndt2['Pb206']/nt['U238'] -
            2*dntdt['Pb206']*dntdt['U238']/nt['U238']^2 -
            out$Pb206U238*d2ndt2['U238']/nt['U238'] +
            2*out$Pb206U238*(dntdt['U238']/nt['U238'])^2
        out$dPb207U235dt <- dntdt['Pb207']/nt['U235'] -
            out$Pb207U235*dntdt['U235']/nt['U235']
        out$d2Pb207U235dt2 <- d2ndt2['Pb207']/nt['U235'] -
            2*dntdt['Pb207']*dntdt['U235']/nt['U235']^2 -
            out$Pb207U235*d2ndt2['U235']/nt['U235'] +
            2*out$Pb207U235*(dntdt['U235']/nt['U235'])^2
        out$dPb206U238dl8 <- 0 # TODO
        out$dPb207U235dl5 <- 0 # TODO
    }
    out
}
forward <- function(tt,d=diseq(),derivative=0){
    if (derivative==0){
        out <- as.vector(d$Q %*% diag(exp(-d$L*tt)) %*% d$Qinv %*% d$atoms)
    } else if (derivative==1){
        out <- as.vector(d$Q %*% diag(exp(-d$L)) %*%
                         diag(exp(-d$L*tt)) %*% d$Qinv %*% d$atoms)
    } else if (derivative==2){
        out <- as.vector(d$Q %*% diag(exp(-d$L)) %*% diag(exp(-d$L)) %*%
                         diag(exp(-d$L*tt)) %*% d$Qinv %*% d$atoms)
    }
    names(out) <- names(d$atoms)
    out
}
reverse <- function(tt,d=diseq(),derivative=FALSE){
    if (derivative==0){
        out <- as.vector(d$Q %*% diag(exp(d$L*tt)) %*% d$Qinv %*% d$atoms)
    } else if (derivative==1){
        out <- as.vector(d$Q %*% diag(exp(d$L)) %*%
                         diag(exp(d$L*tt)) %*% d$Qinv %*% d$atoms)
    } else if (derivative==2){
        out <- as.vector(d$Q %*% diag(exp(d$L)) %*% diag(exp(d$L)) %*%
                         diag(exp(d$L*tt)) %*% d$Qinv %*% d$atoms)        
    }
    names(out) <- names(d$atoms)
    out
}

# from Wendt & Carl (1985, EPSL):
wendt <- function(tt,d=diseq()){
    dd <- geomean(d)
    fact <- 20
    err <- FALSE
    l4 <- settings('lambda','U234')[1]*1000
    l0 <- settings('lambda','Th230')[1]*1000
    l6 <- settings('lambda','Ra226')[1]*1000
    l1 <- settings('lambda','Pa231')[1]*1000    
    # calculate A0, B0, C0 and D0 (Wendt and Carl, 1985)
    # 1. A0
    if (dd$U48$option==0){
        dd$U48$A0 <- 0
    } else if (tt>(fact/l4)){ # expired
        if (dd$U48$option==2) err <- TRUE
        dd$U48$A0 <- dd$U48$x-1
    } else if (dd$U48$option==2){ # back-calculate
        dd$U48$A0 <- max(-1, (dd$U48$x-1)*exp(l4*tt))
    } else { # option 1
        dd$U48$A0 <- dd$U48$x-1
    }
    # 2. B0
    if (dd$ThU$option==0){
        dd$ThU$B0 <- 0
    } else if (tt>(fact/l0)){ # expired
        if (dd$ThU$option==2) err <- TRUE
        dd$ThU$B0 <- dd$ThU$x-1
    } else if (dd$ThU$option==2){
        ThUx <- (dd$ThU$x-1)*exp(l0*tt)
        ThUu <- dd$U48$A0*(1-exp((l0-l4)*tt))*l0/(l0-l4)
        dd$ThU$B0 <- max(-1, ThUu + ThUx)
    } else { # option 1
        dd$ThU$B0 <- dd$ThU$x-1
    }
    # 3. C0
    if (dd$RaU$option==0){
        dd$RaU$C0 <- 0
    } else if (tt>(fact/l6)){ # expired
        if (dd$RaU$option==2) err <- TRUE
        dd$RaU$C0 <- dd$RaU$x-1
    } else if (dd$RaU$option==2){
        RaUx <- (dd$RaU$x-1)*exp(l6*tt)
        RaUth <- dd$ThU$B0*(1-exp((l6-l0)*tt))*l6/(l6-l0)
        RaUu <- dd$U48$A0*(1-exp((l6-l0-l4)*tt))*l6/(l6-l4)
        dd$RaU$C0 <- max(-1, RaUu + RaUth + RaUx)
    } else {
        dd$RaU$C0 <- dd$RaU$x-1
    }
    # 4. D0
    if (dd$PaU$option==0){
        dd$PaU$D0 <- 0
    } else if (tt>(fact/l1)){
        if (dd$PaU$option==2) err <- TRUE
        dd$PaU$D0 <- dd$PaU$x-1
    } else if (dd$PaU$option==2){
        dd$PaU$D0 <- max(-1, (dd$PaU$x-1)*exp(l1*tt))
    } else {
        dd$PaU$D0 <- dd$PaU$x-1
    }
    if (err){
        warning('The measured degree of disequilibrium is ',
                'impossible for a sample of this age.')
    }
    out <- list()
    out$d1 <- d1(tt,dd=dd)
    out$d2 <- d2(tt,dd=dd)
    out$dd1dt <- dd1dt(tt,dd=dd)
    out$dd2dt <- dd2dt(tt,dd=dd)
    out$d2d1dt2 <- d2d1dt2(tt,dd=dd)
    out$d2d2dt2 <- d2d2dt2(tt,dd=dd)
    out$dd1dl5 <- dd1dl5(tt,dd=dd)
    out$dd2dl8 <- dd2dl8(tt,dd=dd)
    out
}
d1 <- function(tt,dd=diseq()){
    l1 <- settings('lambda','Pa231')[1]*1000
    l5 <- settings('lambda','U235')[1]
    dd$PaU$D0*(l5/l1)*(1-exp(-l1*tt))
}
d2 <- function(tt,dd=diseq()){
    K <- K1234(dd)
    K$K1*exp(-K$l4*tt) + K$K2*exp(-K$l0*tt) + K$K3*exp(-K$l6*tt) + K$K4
}
K1234 <- function(dd=diseq()){
    l8 <- settings('lambda','U238')[1]
    l4 <- settings('lambda','U234')[1]*1000
    l0 <- settings('lambda','Th230')[1]*1000
    l6 <- settings('lambda','Ra226')[1]*1000
    out <- list(l8=l8,l4=l4,l0=l0,l6=l6)
    out$K1 <- -dd$U48$A0*l8*l0*l6/(l4*(l0-l4)*(l6-l4))
    out$K2 <- (l8*l6/(l6-l0))*(dd$U48$A0/(l0-l4)-dd$ThU$B0/l0)
    out$K3 <- (l8/(l6-l0))*(dd$ThU$B0-l0*dd$U48$A0/(l6-l4))-dd$RaU$C0*l8/l6
    out$K4 <- dd$U48$A0*l8/l4 + dd$ThU$B0*l8/l0 + dd$RaU$C0*l8/l6
    out
}
dd1dl5 <- function(tt,dd=diseq()){
    l5 <- settings('lambda','U235')[1]
    d1(tt,dd=dd)/l5
}
dd2dl8 <- function(tt,dd=diseq()){
    l8 <- settings('lambda','U238')[1]
    d2(tt,dd=dd)/l8
}
dd1dt <- function(tt,dd=diseq()){
    l5 <- settings('lambda','U235')[1]
    l1 <- settings('lambda','Pa231')[1]*1000
    l5*dd$PaU$D0*exp(-l1*tt)
}
d2d1dt2 <- function(tt,dd=diseq()){
    l1 <- settings('lambda','Pa231')[1]*1000
    -l1*dd1dt(tt,dd=dd)
}
dd2dt <- function(tt,dd=diseq()){
    K <- K1234(dd)
    out <- - K$l4*K$K1*exp(-K$l4*tt) - K$l0*K$K2*exp(-K$l0*tt) -
        K$l6*K$K3*exp(-K$l6*tt)
    out
}
d2d2dt2 <- function(tt,dd=diseq()){
    K <- K1234(dd)
    out <- K$K1*exp(-K$l4*tt)*K$l4^2 + K$K2*exp(-K$l0*tt)*K$l0^2 +
        K$K3*exp(-K$l6*tt)*K$l6^2
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
    D <- wendt(tt=t.76,d=d)
    r75 <- exp(l5*t.76)-1+D$d1
    r68 <- exp(l8*t.76)-1+D$d2
    dr75dt <- l5*exp(l5*t.76)+D$dd1dt
    dr68dt <- l8*exp(l8*t.76)+D$dd2dt
    2*(r75/(U*r68)-x)*(dr75dt*r68-r75*dr68dt)/(U*r68^2)
}
dmf76dl5 <- function(x,t.76,d=diseq()){
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    U <- iratio('U238U235')[1]
    D <- wendt(tt=t.76,d=d)
    r75 <- exp(l5*t.76)-1+D$d1
    r68 <- exp(l8*t.76)-1+D$d2
    dr75dl5 <- t.76*exp(l5*t.76)+D$dd1dl5
    2*(r75/(U*r68)-x)/(U*r68)
}
dmf76dl8 <- function(x,t.76,d=diseq()){
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    U <- iratio('U238U235')[1]
    D <- wendt(tt=t.76,d=d)
    r75 <- exp(l5*t.76)-1+D$d1
    r68 <- exp(l8*t.76)-1+D$d2
    dr68dl8 <- l8*exp(l8*t.76)+D$dd2dt
    2*(r75/(U*r68)-x)*(-r75*dr68dl8)/(U*r68^2)
}
dmf76dU <- function(x,t.76,d=diseq()){
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    U <- iratio('U238U235')[1]
    D <- wendt(tt=t.76,d=d)
    r75 <- exp(l5*t.76)-1+D$d1
    r68 <- exp(l8*t.76)-1+D$d2
    -2*(r75/(U*r68)-x)*r75/(r68*U^2)
}
