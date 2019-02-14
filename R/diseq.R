diseq <- function(U48=1,Th0U8=1,Ra6U8=1,Pa1U5=1){
    out <- list(U48=U48,Th0U8=Th0U8,Ra6U8=Ra6U8,Pa1U5=Pa1U5,corr=FALSE)
    if (U48!=1 | Th0U8!=1 | Ra6U8!=1 | Pa1U5!=1) out$corr <- TRUE
    class(out) <- 'diseq'
    out
}

# from Wendt & Carl (1985, EPSL):
wendt <- function(tt,d=diseq()){
    out <- list(d1=0,d2=0,dd1dt=0,dd2dt=0,
                d2d1dt2=0,d2d2dt2=0,dd1dl5=0,dd2dl8=0)
    if (d$corr){
        out$d1 <- d1(tt,Pa1U5=d$Pa1U5)
        out$d2 <- d2(tt,U48=d$U48,Th0U8=d$Th0U8,Ra6U8=d$Ra6U8)
        out$dd1dt <- dd1dt(tt,Pa1U5=d$Pa1U5)
        out$dd2dt <- dd2dt(tt,U48=d$U48,Th0U8=d$Th0U8,Ra6U8=d$Ra6U8)
        out$d2d1dt2 <- d2d1dt2(tt,Pa1U5=d$Pa1U5)
        out$d2d2dt2 <- d2d2dt2(tt,U48=d$U48,Th0U8=d$Th0U8,Ra6U8=d$Ra6U8)
        out$dd1dl5 <- dd1dl5(tt,Pa1U5=d$Pa1U5)
        out$dd2dl8 <- dd2dl8(tt,U48=d$U48,Th0U8=d$Th0U8,Ra6U8=d$Ra6U8)
    }
    out
}
d1 <- function(tt,Pa1U5=0){
    l1 <- settings('lambda','Pa231')[1]*1000
    l5 <- settings('lambda','U235')[1]
    D0 <- Pa1U5 - 1
    D0*(l5/l1)*exp(l5*tt)*(1-exp(-l1*tt))
}
d2 <- function(tt,U48=1,Th0U8=0,Ra6U8=0){
    A0 <- U48 - 1
    B0 <- Th0U8 - 1
    C0 <- Ra6U8 - 1
    l6 <- settings('lambda','Ra226')[1]*1000
    l0 <- settings('lambda','Th230')[1]*1000
    l4 <- settings('lambda','U234')[1]*1000
    l8 <- settings('lambda','U238')[1]
    K1 <- -A0*l8*l0*l6/(l4*(l0-l4)*(l6-l4))
    K2 <- (l8*l6/(l6-l0))*(A0/(l0-l4)-B0/l0)
    K3 <- (l8/(l6-l0))*(B0-l0*A0/(l6-l4))-C0*l8/l6
    K4 <- A0*l8/l4 + B0*l8/l0 + C0*l8/l6
    exp(l8*tt)*(K1*exp(-l4*tt)+K2*exp(-l0*tt)+K3*exp(-l6*tt)+K4)
}
dd1dl5 <- function(tt,Pa1U5=0){
    l5 <- settings('lambda','U235')[1]
    d1(tt,Pa1U5=Pa1U5)*(tt+1/l5)
}
dd2dl8 <- function(tt,U48=1,Th0U8=0,Ra6U8=0){
    l8 <- settings('lambda','U238')[1]
    d2(tt,U48=U48,Th0U8=Th0U8,Ra6U8=Ra6U8)*(tt+1/l8)
}
dd1dt <- function(tt,Pa1U5=0){
    D0 <- Pa1U5 - 1
    l5 <- settings('lambda','U235')[1]
    l1 <- settings('lambda','Pa231')[1]*1000
    l5*d1(tt,Pa1U5=Pa1U5) + D0*l5*exp((l5-l1)*tt)
}
d2d1dt2 <- function(tt,Pa1U5=0){
    D0 <- Pa1U5 - 1
    l5 <- settings('lambda','U235')[1]
    l1 <- settings('lambda','Pa231')[1]*1000
    l5*l5*d1(tt,Pa1U5=Pa1U5) + D0*l5*(l5-l1)*exp((l5-l1)*tt)
}
dd2dt <- function(tt,U48=1,Th0U8=0,Ra6U8=0){
    A0 <- U48 - 1
    B0 <- Th0U8 - 1
    C0 <- Ra6U8 - 1    
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
    out
}
d2d2dt2 <- function(tt,U48=1,Th0U8=0,Ra6U8=0){
    A0 <- U48 - 1
    B0 <- Th0U8 - 1
    C0 <- Ra6U8 - 1    
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
