get.ThU.age <- function(U234U238,sU234U238,Th230U238,sTh230U238,cov4808,
                        exterr=TRUE,cor=FALSE){
    l0 <- lambda('Th230')
    l4 <- lambda('U234')
    a <- U234U238
    A <- Th230U238
    sa <- sU234U238
    sA <- sTh230U238
    covAa <- cov4808
    fit <- stats::optim(0,fn=ThU.misfit,gr=ThU.gr,method='BFGS',
                        A=A,a=a,l0=l0,l4=l4)
    tt <- fit$par
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
    ns <- length(x)
    d <- data2evolution(x)
    if (i2i){
        osmond <- data2tit.ThU(x,osmond=TRUE)
        fit <- titterington(osmond)
    }
    out <- matrix(0,ns,5)
    colnames(out) <- c('t','s[t]','48_0','s[48_0]','cov[t,48_0]')
    if (cor) colnames(out)[5] <- 'cor[t,48_0]'
    for (j in 1:ns){
        a <- d[j,'U234U238']
        sa <- d[j,'errU234U238']
        A <- d[j,'Th230U238']
        sA <- d[j,'errTh230U238']
        covAa <- d[j,'cov']
        if (i2i){ # project onto 08-48 plane
            a0 <- osmond[j,'Y'] - fit$par['b']*osmond[j,'X']
            A0 <- osmond[j,'Z'] - fit$par['B']*osmond[j,'X']
            out[j,] <- get.ThU.age(a0,sa,A0,sA,covAa,exterr=exterr,cor=cor)
        } else {
            out[j,] <- get.ThU.age(a,sa,A,sA,covAa,exterr=exterr,cor=cor)
        }
    }
    if (!is.na(sigdig)){
        out[,c(1,2)] <- roundit(out[,1],out[,2],sigdig=sigdig)
        out[,c(3,4)] <- roundit(out[,3],out[,4],sigdig=sigdig)
        out[,5] <- signif(out[,5],sigdig)
    }
    if (!is.na(i)) out <- out[i,]
    out
}
