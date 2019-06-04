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
#' @return
#' a list with the following items: \code{option} and (\code{U48},
#' \code{Th08}, \code{Ra6U8}, \code{Pa1U8}) [if \code{option=1} or
#' \code{option=2}] and (\code{fThU}, \code{RaU}, \code{PaU}) [if
#' \code{option=3}].
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
diseq <- function(option=0,
                  U48=1,Th0U8=1,Ra6U8=1,Pa1U5=1,
                  fThU=1,fRaU=1,fPaU=1){
    out <- list(option=option)
    if (option%in%c(1,2)){
        out$U48=U48
        out$Th0U8=Th0U8
        out$Ra6U8=Ra6U8
        out$Pa1U5=Pa1U5        
    } else if (option==3){
        out$fThU = fThU
        out$fRaU = fRaU
        out$fPaU = fPaU
    }
    out
}

# from Wendt & Carl (1985, EPSL):
wendt <- function(tt,d=diseq()){
    dd <- d
    if (d$option==2){
        l4 <- settings('lambda','U234')[1]*1000
        l0 <- settings('lambda','Th230')[1]*1000
        l6 <- settings('lambda','Ra226')[1]*1000
        l1 <- settings('lambda','Pa231')[1]*1000
        # the next 6 lines are safety measures against bad input
        fact <- 20
        if (tt>(fact/l6) & d$Ra6U8!=1) ttt <- fact/l6
        else if (tt>(fact/l1) & d$Pa1U5!=1) ttt <- fact/l1
        else if (tt>(fact/l0) & d$Th0U8!=1) ttt <- fact/l0
        else if (tt>(fact/l4) & d$U48!=1) ttt <- fact/l4
        else ttt <- tt
        if (d$U48!=1) dd$U48 <- max(0, 1 + (d$U48-1)*exp(l4*ttt))
        if (d$Th0U8!=1) dd$Th0U8 <- max(0, 1 + (d$Th0U8-1)*exp(l0*ttt))
        if (d$Ra6U8!=1) dd$Ra6U8 <- max(0, 1 + (d$Ra6U8-1)*exp(l6*ttt))
        if (d$Pa1U5!=1) dd$Pa1U5 <- max(0, 1 + (d$Pa1U5-1)*exp(l1*ttt))
    }
    out <- list(d1=0,d2=0,dd1dt=0,dd2dt=0,
                d2d1dt2=0,d2d2dt2=0,dd1dl5=0,dd2dl8=0)
    if (d$option>0){
        out$d1 <- d1(tt,dd=dd)
        out$d2 <- d2(tt,dd=dd)
        out$dd1dt <- dd1dt(tt,dd=dd)
        out$dd2dt <- dd2dt(tt,dd=dd)
        out$d2d1dt2 <- d2d1dt2(tt,dd=dd)
        out$d2d2dt2 <- d2d2dt2(tt,dd=dd)
        out$dd1dl5 <- dd1dl5(tt,dd=dd)
        out$dd2dl8 <- dd2dl8(tt,dd=dd)
    }
    out
}
d1 <- function(tt,dd=diseq()){
    l1 <- settings('lambda','Pa231')[1]*1000
    l5 <- settings('lambda','U235')[1]
    if (dd$option<3){
        D0 <- dd$Pa1U5 - 1
        out <- D0*(l5/l1)*exp(l5*tt)*(1-exp(-l1*tt))
    } else {
        out <- (l5/l1)*(dd$fPaU-1)
    }
    out
}
d2 <- function(tt,dd=diseq()){
    l4 <- settings('lambda','U234')[1]*1000
    l6 <- settings('lambda','Ra226')[1]*1000
    l0 <- settings('lambda','Th230')[1]*1000
    l8 <- settings('lambda','U238')[1]
    if (dd$option<3){
        A0 <- dd$U48 - 1
        B0 <- dd$Th0U8 - 1
        C0 <- dd$Ra6U8 - 1
        K1 <- -A0*l8*l0*l6/(l4*(l0-l4)*(l6-l4))
        K2 <- (l8*l6/(l6-l0))*(A0/(l0-l4)-B0/l0)
        K3 <- (l8/(l6-l0))*(B0-l0*A0/(l6-l4))-C0*l8/l6
        K4 <- A0*l8/l4 + B0*l8/l0 + C0*l8/l6
        out <- exp(l8*tt)*(K1*exp(-l4*tt)+K2*exp(-l0*tt)+K3*exp(-l6*tt)+K4)
    } else {
        out <- (l8/l0)*(dd$fThU-1) + (l8/l6)*(dd$fRaU-1)
    }
    out
}
dd1dl5 <- function(tt,dd=diseq()){
    if (dd$option<3){
        l5 <- settings('lambda','U235')[1]
        out <- d1(tt,dd=dd)*(tt+1/l5)
    } else {
        l1 <- settings('lambda','Pa231')[1]*1000
        out <- (dd$fPaU-1)/l1
    }
    out
}
dd2dl8 <- function(tt,dd=diseq()){
    if (dd$option<3){
        l8 <- settings('lambda','U238')[1]
        out <- d2(tt,dd=dd)*(tt+1/l8)
    } else {
        l0 <- settings('lambda','Th230')[1]*1000
        l6 <- settings('lambda','Ra226')[1]*1000
        out <- (dd$fThU-1)/l0 + (dd$fRaU-1)/l6
    }
    out
}
dd1dt <- function(tt,dd=diseq()){
    if (dd$option<3){
        l5 <- settings('lambda','U235')[1]
        l1 <- settings('lambda','Pa231')[1]*1000
        D0 <- dd$Pa1U5 - 1
        out <- l5*d1(tt,dd=dd) + D0*l5*exp((l5-l1)*tt)
    } else {
        out <- 0
    }
    out
}
d2d1dt2 <- function(tt,dd=diseq()){
    if (dd$option<3){
        D0 <- dd$Pa1U5 - 1
        l5 <- settings('lambda','U235')[1]
        l1 <- settings('lambda','Pa231')[1]*1000
        out <- l5*l5*d1(tt,dd=dd) + D0*l5*(l5-l1)*exp((l5-l1)*tt)
    } else {
        out <- 0
    }
    out
}
dd2dt <- function(tt,dd=diseq()){
    if (dd$option<3){
        A0 <- dd$U48 - 1
        B0 <- dd$Th0U8 - 1
        C0 <- dd$Ra6U8 - 1    
        l6 <- settings('lambda','Ra226')[1]*1000
        l0 <- settings('lambda','Th230')[1]*1000
        l4 <- settings('lambda','U234')[1]*1000
        l8 <- settings('lambda','U238')[1]
        K1 <- -A0*l8*l0*l6/(l4*(l0-l4)*(l6-l4))
        K2 <- (l8*l6/(l6-l0))*(A0/(l0-l4)-B0/l0)
        K3 <- (l8/(l6-l0))*(B0-l0*A0/(l6-l4))-C0*l8/l6
        K4 <- A0*l8/l4 + B0*l8/l0 + C0*l8/l6
        out <- (l8-l4)*K1*exp((l8-l4)*tt) + (l8-l0)*K2*exp((l8-l0)*tt) +
            (l8-l6)*K3*exp((l8-l6)*tt) + l8*K4*exp(l8*tt)
    } else {
        out <- 0
    }
    out
}
d2d2dt2 <- function(tt,dd=diseq()){
    if (dd$option<3){
        A0 <- dd$U48 - 1
        B0 <- dd$Th0U8 - 1
        C0 <- dd$Ra6U8 - 1
        l6 <- settings('lambda','Ra226')[1]*1000
        l0 <- settings('lambda','Th230')[1]*1000
        l4 <- settings('lambda','U234')[1]*1000
        l8 <- settings('lambda','U238')[1]
        K1 <- -A0*l8*l0*l6/(l4*(l0-l4)*(l6-l4))
        K2 <- (l8*l6/(l6-l0))*(A0/(l0-l4)-B0/l0)
        K3 <- (l8/(l6-l0))*(B0-l0*A0/(l6-l4))-C0*l8/l6
        K4 <- A0*l8/l4 + B0*l8/l0 + C0*l8/l6
        out <- (l8-l4)*(l8-l4)*K1*exp((l8-l4)*tt) +
            (l8-l0)*(l8-l0)*K2*exp((l8-l0)*tt) +
            (l8-l6)*(l8-l6)*K3*exp((l8-l6)*tt) + K4*l8*exp(l8*tt)
    } else {
        out <- 0
    }
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
