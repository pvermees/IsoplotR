get.ThU.age <- function(U234U238,sU234U238,Th230U238,sTh230U238,cov4808,exterr=TRUE){
    a <- U234U238
    sa <- sU234U238
    A <- Th230U238
    sA <- sTh230U238
    covAa <- cov4808
    l4 <- lambda('U234')
    l0 <- lambda('Th230')
    fit <- stat::optim(0.1,fn=ThU.misfit,gr=ThU.gr,method='BFGS',A=A,a=a,l0=l0,l4=l4)
    tt <- fit$par
    D <- l0[1]*(exp(-l0[1]*tt)+(a-1)*exp((l4[1]-l0[1])*tt))
    k1 <- k1(tt,l0,l4)
    st <- sqrt((sA^22 + (k1*sa)^2 + 2*k1*covAa)/D^2)
    1000*c(tt,st)
}

k1 <- function(tt,l0,l4){
    (1-exp((l4[1]-l0[1])*tt))*l0[1]/(l4[1]-l0[1])
}
ThU.misfit <- function(tt,A,a,l0,l4){
    (1-exp(-l0[1]*tt)-(a-1)*k1(tt,l0,l4) - A)^2
}

ThU.gr <- function(tt,A,a,l0,l4){
    dk1.dtt <- l0[1]*exp((l4[1]-l0[1])*tt)
    dmisfit.dtt <- 2*(1-exp(-l0[1]*tt)-(a-1)*k1(tt,l0,l4) - A) *
        (l0[1]*exp(-l0[1]*tt) - (a-1)*dk1.dtt)
    dmisfit.dtt
}


