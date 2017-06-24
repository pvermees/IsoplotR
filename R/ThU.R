get.ThU.age <- function(U234U238,sU234U238,Th230U238,sTh230U238,cov4808,exterr=TRUE){
    a <- U234U238
    sa <- sU234U238
    A <- Th230U238
    sA <- sTh230U238
    covAa <- cov4808
    l4 <- lambda('U234')
    l0 <- lambda('Th230')
    fit <- stats::optim(0,fn=ThU.misfit,gr=ThU.gr,method='BFGS',A=A,a=a,l0=l0,l4=l4)
    tt <- fit$par
    l40 <- l4[1]-l0[1]
    el40t <- exp(l40*tt)
    dk1.dl0 <- l4[1]*(1-el40t)/l40^2 + l0[1]*tt*el40t/l40
    dk1.dl4 <- l0[1]*(el40t-1)/l40^2 - l0[1]*tt*el40t/l40
    dk1.dt <- -l0[1]*el40t
    dD.da <- k1(tt,l0,l4)
    dD.dA <- 1
    dD.dl0 <- (a-1)*dk1.dl0 - tt*exp(-l0[1]*tt)
    dD.dl4 <- (a-1)*dk1.dl4
    dD.dt <- (a-1)*dk1.dt - l0[1]*exp(-l0[1]*tt)
    J <- matrix(0,1,4)
    J[1,1] <- -dD.da/dD.dt
    J[1,2] <- -dD.dA/dD.dt
    if (exterr){
        J[1,3] <- dD.dl0/dD.dt
        J[1,4] <- -dD.dl4/dD.dt
    }
    E <- matrix(0,4,4)
    diag(E) <- c(sa,sA,l0[2],l4[2])^2
    E[1,2] <- covAa
    E[2,1] <- E[1,2]
    st <- sqrt(J %*% E %*% t(J))
    # The above gives the same result as equation 6 of 
    # Ludwig and Titterington (2004), with a fixed typo:
    # D <- l0[1]*(exp(-l0[1]*tt)+(a-1)*exp((l4[1]-l0[1])*tt))
    # k1 <- k1(tt,l0,l4)
    # st <- sqrt((sA^22 + (k1*sa)^2 - 2*k1*covAa)/D^2)
    c(tt,st)
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


