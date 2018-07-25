get.ThU.age <- function(Th230U238,sTh230U238,
                        U234U238=1,sU234U238=0,
                        cov4808=0,exterr=TRUE,cor=FALSE){
    l0 <- lambda('Th230')
    l4 <- lambda('U234')
    a <- U234U238
    A <- Th230U238
    sa <- sU234U238
    sA <- sTh230U238
    covAa <- cov4808
    fit <- fit <- stats::optimize(ThU.misfit,interval=c(0,2000),A=A,a=a,l0=l0,l4=l4)
    tt <- fit$minimum
    a0 <- 1+(a-1)*exp(l4[1]*tt)
    l40 <- l4[1]-l0[1]
    el40t <- exp(l40*tt)
    el4t <- exp(l4[1]*tt)
    dk1.dl0 <- l4[1]*(1-el40t)/l40^2 + l0[1]*tt*el40t/l40
    dk1.dl4 <- l0[1]*(el40t-1)/l40^2 - l0[1]*tt*el40t/l40
    dk1.dt <- -l0[1]*el40t
    dD.da <- k1(tt,l0,l4)
    dD.dA <- 1
    dD.dl0 <- (a-1)*dk1.dl0 - tt*exp(-l0[1]*tt)
    dD.dl4 <- (a-1)*dk1.dl4
    dD.dt <- (a-1)*dk1.dt - l0[1]*exp(-l0[1]*tt)
    dt.da <- -dD.da/dD.dt
    dt.dA <- -dD.dA/dD.dt
    dt.dl0 <- -dD.dl0/dD.dt
    dt.dl4 <- -dD.dl4/dD.dt
    da0.da <- el4t + (a-1)*l4[1]*dt.da*el4t
    da0.dA <- (a-1)*l4[1]*dt.dA*el4t
    da0.dl0 <- (a-1)*l4[1]*dt.dl0*el4t
    da0.dl4 <- (a-1)*tt*el4t + (a-1)*l4[1]*dt.dl4*el4t
    J <- matrix(0,2,4)
    J[1,1] <- dt.da
    J[1,2] <- dt.dA
    J[2,1] <- da0.da
    J[2,2] <- da0.dA
    if (exterr){
        J[1,3] <- dt.dl0
        J[1,4] <- dt.dl4
        J[2,3] <- da0.dl0
        J[2,4] <- da0.dl4
    }
    E <- matrix(0,4,4)
    diag(E) <- c(sa,sA,l0[2],l4[2])^2
    E[1,2] <- covAa
    E[2,1] <- E[1,2]
    covmat <- J %*% E %*% t(J)
    st <- sqrt(covmat[1,1])
    sa0 <- sqrt(covmat[2,2])
    out <- c(tt,st,a0,sa0,covmat[1,2])
    names(out) <- c('t','s[t]','48_0','s[48_0]','cov[t,48_0]')
    if (cor){
        out[5] <- covmat[1,2]*st*sa0
        names(out)[5] <- 'cor[t,48_0]'
    }
    out
}

get.Th230Th232_0x <- function(tt,Th230Th232,errTh230Th232=0){
    l0 <- settings('lambda','Th230')[1]
    Th230Th232_0x <- Th230Th232 * exp(l0*tt)
    errTh230Th232_0x <- errTh230Th232/Th230Th232
    c(Th230Th232_0x,errTh230Th232_0x)
}

get.Th230Th232 <- function(tt,Th230Th232_0x,U238Th232){
    l0 <- settings('lambda','Th230')[1]
    Th230Th232_0x * exp(-l0*tt) + U238Th232 * (1 - exp(-l0*tt))
}

get.U238Th232 <- function(tt,Th230Th232_0x,Th230Th232){
    l0 <- settings('lambda','Th230')[1]
    (Th230Th232 - Th230Th232_0x * exp(-l0*tt)) / (1 - exp(-l0*tt))
}

k1 <- function(tt,l0,l4){
    (1-exp((l4[1]-l0[1])*tt))*l0[1]/(l4[1]-l0[1])
}
ThU.misfit <- function(tt,A,a,l0,l4){
    (1-exp(-l0[1]*tt)-(a-1)*k1(tt,l0,l4) - A)^2
}
# not used anymore but kept for possible use in the future
ThU.gr <- function(tt,A,a,l0,l4){
    dk1.dtt <- l0[1]*exp((l4[1]-l0[1])*tt)
    dmisfit.dtt <- 2*(1-exp(-l0[1]*tt)-(a-1)*k1(tt,l0,l4) - A) *
        (l0[1]*exp(-l0[1]*tt) - (a-1)*dk1.dtt)
    dmisfit.dtt
}

# converts format 1 to format 2 and vice versa
ThU.convert <- function(x){
    labels <- c('X','sX','Y','sY','Z','sZ','rXY','rXZ','rYZ')
    out <- rep(0,9)
    names(x) <- labels
    names(out) <- labels
    out['X'] <- 1/x['X']
    out['Y'] <- x['Y']/x['X']
    out['Z'] <- x['Z']/x['X']
    J <- matrix(0,3,3)
    J[1,1] <- -1/x['X']^2
    J[1,2] <- 0
    J[1,3] <- 0
    J[2,1] <- -x['Y']/x['X']^2
    J[2,2] <- 1/x['X']
    J[2,3] <- 0
    J[3,1] <- -x['Z']/x['X']^2
    J[2,3] <- 0
    J[3,3] <- 1/x['X']
    E <- cor2cov3(x['sX'],x['sY'],x['sZ'],x['rXY'],x['rXZ'],x['rYZ'])
    covmat <- J %*% E %*% t(J)
    cormat <- stats::cov2cor(covmat)
    out['sX'] <- sqrt(covmat[1,1])
    out['sY'] <- sqrt(covmat[2,2])
    out['sZ'] <- sqrt(covmat[3,3])
    out['rXY'] <- cormat[1,2]
    out['rXZ'] <- cormat[1,3]
    out['rYZ'] <- cormat[2,3]
    out
}

get.Th230U238 <- function(tt,U234U238_0){
    l4 <- lambda('U234')
    l0 <- lambda('Th230')
    a <- get.U234U238(tt,U234U238_0)
    1 - exp(-l0[1]*tt) - (a-1)*k1(tt,l0,l4)
}

get.U234U238 <- function(tt,U234U238_0){
    l4 <- lambda('U234')
    1 + (U234U238_0-1)*exp(-l4[1]*tt)
}

ThU.age <- function(x,exterr=FALSE,i=NA,i2i=FALSE,sigdig=NA,cor=TRUE){
    if (x$format %in% c(1,2))
        out <- get.ThU.age.corals(x,exterr=exterr,i=i,i2i=i2i,sigdig=sigdig,cor=cor)
    else
        out <- get.ThU.age.volcanics(x,exterr=exterr,i=i,i2i=i2i,sigdig=sigdig)
    out
}

get.ThU.age.corals <- function(x,exterr=FALSE,i=NA,i2i=FALSE,sigdig=NA,cor=TRUE){
    ns <- length(x)
    if (ns<2) i2i <- FALSE
    osmond <- data2tit.ThU(x,osmond=TRUE)
    if (i2i) fit <- titterington(osmond)
    out <- matrix(0,ns,5)
    colnames(out) <- c('t','s[t]','48_0','s[48_0]','cov[t,48_0]')
    if (cor) colnames(out)[5] <- 'cor[t,48_0]'
    Th0Th2_d<- settings('iratio','Th230Th232')[1]
    sTh0Th2_d <- settings('iratio','Th230Th232')[2]
    U4U8_d <- settings('iratio','U234U238')[1]
    sU4U8_d <- settings('iratio','U234U238')[2]
    cov4808 <- data2evolution(x)[,'cov']
    for (j in 1:ns){
        Th2U8 <- osmond[j,'X']
        sTh2U8 <- osmond[j,'sX']
        U4U8 <- osmond[j,'Y']
        sU4U8 <- osmond[j,'sY']
        Th0U8 <- osmond[j,'Z']
        sTh0U8 <- osmond[j,'sZ']
        if (i2i){ # project onto 08-48 plane
            U4Th2 <- fit$par['b']
            U4U8_d <- U4Th2 * Th2U8
            sU4U8_d <- 0 # don't propagate projection error
            Th0Th2 <- fit$par['B']
            Th0U8_d <- Th0Th2 * Th2U8
            sTh0U8_d <- 0 # don't propagate projection error
        } else {
            Th0U8_d <- Th0Th2_d * Th2U8
            if (Th0U8_d > 0)
                sTh0U8_d <- Th0U8_d * sqrt( (sTh0Th2_d/Th0Th2_d)^2 + (sTh2U8/Th2U8)^2 )
            else
                sTh0U8_d <- 0
        }
        U4U8_0 <- U4U8 - U4U8_d
        if (U4U8_d > 0)
            sU4U8_0 <- U4U8_0 * sqrt( (sU4U8/U4U8)^2 + (sU4U8_d/U4U8_d)^2 )
        else
            sU4U8_0 <- sU4U8
        Th0U8_0 <- Th0U8 - Th0U8_d
        if (Th0U8_d > 0)
            sTh0U8_0 <- Th0U8_0 * sqrt( (sTh0U8/Th0U8)^2 + (sTh0U8_d/Th0U8_d)^2 )
        else
            sTh0U8_0 <- sTh0U8
        out[j,] <- get.ThU.age(Th0U8_0,sTh0U8_0,U4U8_0,sU4U8_0,
                               cov4808[j],exterr=exterr,cor=cor)
    }
    if (!is.na(sigdig)){
        temp <- out # store numerical values because roundit creates text
        temp[,c(1,2)] <- roundit(out[,1],out[,2],sigdig=sigdig)
        temp[,c(3,4)] <- roundit(out[,3],out[,4],sigdig=sigdig)
        temp[,5] <- signif(out[,5],sigdig)
        out <- temp
    }
    if (!is.na(i)) out <- out[i,]
    out
}

get.ThU.age.volcanics <- function(x,exterr=FALSE,i=NA,i2i=FALSE,sigdig=NA){
    ns <- length(x)
    d <- data2evolution(x,isochron=i2i)
    Th230U238 <- (d$x[,'Th230Th232'] - d$Th230Th232_0)/d$x[,'U238Th232']
    d08.d02 <- 1/d$x[,'U238Th232']
    d08.d82 <- -(d$x[,'Th230Th232'] - d$Th230Th232_0)/d$x[,'U238Th232']^2
    sTh230U238 <- sqrt( (d08.d02*d$x[,'errTh230Th232'])^2 +
                        (d08.d82*d$x[,'errU238Th232'])^2 )
    out <- matrix(0,ns,2)
    colnames(out) <- c('t','s[t]')
    for (j in 1:ns){
        tst <- get.ThU.age(Th230U238[j],sTh230U238[j],exterr=exterr)
        out[j,] <- tst[c('t','s[t]')]
    }
    if (!is.na(sigdig)) out <- roundit(out[,1],out[,2],sigdig=sigdig)
    if (!is.na(i)) out <- out[i,]
    out
}
