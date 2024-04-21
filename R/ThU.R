ThU_age <- function(x,exterr=FALSE,i=NULL,Th0i=0,
                    cor=TRUE,omit4c=NULL){
    if (x$format %in% c(1,2)){
        out <- get_ThU_age_corals(x,exterr=exterr,i=i,cor=cor,
                                  Th0i=Th0i,omit4c=omit4c)
    } else {
        out <- get_ThU_age_volcanics(x,exterr=exterr,i=i,
                                     Th0i=Th0i,omit4c=omit4c)
    }
    out
}

get_ThU_age_corals <- function(x,exterr=FALSE,i=NULL,
                               cor=TRUE,Th0i=0,omit4c=NULL){
    td <- data2tit.ThU(x,osmond=TRUE,generic=FALSE) # 2/8 - 4/8 - 0/8
    if (Th0i==2){
        td <- Th230correction_measured_detritus(td,Th02U48=x$Th02U48)
    } else {
        if (Th0i==1){
            td <- Th230correction_isochron(td,omit4c=omit4c)
        }
    }
    if (Th0i==3) Th02i <- x$Th02i
    else Th02i <- c(0,0)
    ns <- length(x)
    out <- matrix(0,ns,5)
    for (j in 1:ns){
        out[j,] <- get_ThU_age(Th230U238=td[j,'Th230U238'],sTh230U238=td[j,'sTh230U238'],
                               U234U238=td[j,'U234U238'],sU234U238=td[j,'sU234U238'],
                               cov4808=td[j,'rYZ']*td[j,'sU234U238']*td[j,'sTh230U238'],
                               Th232U238=td[j,'Th232U238'],sTh232U238=td[j,'sTh232U238'],
                               Th230Th232i=Th02i[1],sTh230Th232i=Th02i[2],
                               exterr=exterr,cor=cor)
    }
    if (!is.null(i)) out <- out[i,]
    colnames(out) <- c('t','s[t]','48_0','s[48_0]','cov[t,48_0]')
    out
}

get_ThU_age_volcanics <- function(x,exterr=FALSE,i=NULL,Th0i=0,omit4c=NULL){
    ns <- length(x)
    d <- data2york(x,type=2,generic=FALSE)
    if (Th0i==1){
        fit <- isochron.ThU(x,type=2,plot=FALSE,exterr=FALSE,omit=omit4c)
        Th02 <- fit$b[1]
    } else if (Th0i==2){
        Th02 <- (1-d[,'Th230U238'])/(1/x$U8Th2-d[,'Th232U238'])
    } else {
        Th02 <- 0
    }
    d[,'Th230U238'] <- d[,'Th230U238'] - Th02*d[,'Th232U238']
    out <- matrix(0,ns,2)
    for (j in 1:ns){
        tst <- get_ThU_age(d[j,'Th230U238'],d[j,'errTh230U238'],exterr=exterr)
        out[j,] <- tst[c('t','s[t]')]
    }
    if (!is.null(i)) out <- out[i,]
    colnames(out) <- c('t','s[t]')
    out
}

get_ThU_age <- function(Th230U238,sTh230U238,
                        U234U238=1,sU234U238=0,cov4808=0,
                        Th232U238=0,sTh232U238=0,
                        Th230Th232i=0,sTh230Th232i=0,
                        exterr=FALSE,cor=FALSE,jacobian=FALSE){
    l0 <- lambda('Th230')
    l4 <- lambda('U234')
    a <- U234U238
    A <- Th230U238
    sa <- sU234U238
    sA <- sTh230U238
    covAa <- cov4808
    fit <- stats::optimize(ThU_misfit,interval=c(0,2000),
                           A=A,a=a,l0=l0,l4=l4,
                           Th2U8=Th232U238,Th02i=Th230Th232i)
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
    dD.dTh2U8 <- -exp(-l0[1]*tt)*Th230Th232i
    dD.dTh02i <- -exp(-l0[1]*tt)*Th232U238
    dt.da <- -dD.da/dD.dt
    dt.dA <- -dD.dA/dD.dt
    dt.dl0 <- -dD.dl0/dD.dt
    dt.dl4 <- -dD.dl4/dD.dt
    dt.dTh2U8 <- -dD.dTh2U8/dD.dt
    dt.dTh02i <- -dD.dTh02i/dD.dt
    da0.da <- el4t + (a-1)*l4[1]*dt.da*el4t
    da0.dA <- (a-1)*l4[1]*dt.dA*el4t
    da0.dl0 <- (a-1)*l4[1]*dt.dl0*el4t
    da0.dl4 <- (a-1)*tt*el4t + (a-1)*l4[1]*dt.dl4*el4t
    J <- matrix(0,2,6)
    J[1,1] <- dt.da
    J[1,2] <- dt.dA
    J[2,1] <- da0.da
    J[2,2] <- da0.dA
    if (exterr){
        J[1,3] <- dt.dl0
        J[1,4] <- dt.dl4
        J[1,5] <- dt.dTh2U8
        J[1,6] <- dt.dTh02i
        J[2,3] <- da0.dl0
        J[2,4] <- da0.dl4
    }
    E <- matrix(0,6,6)
    diag(E) <- c(sa,sA,l0[2],l4[2],sTh232U238,sTh230Th232i)^2
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
    if (jacobian){
        out <- c(out,J[1,1:4])
        names(out)[6:9] <- c('dt.d48','dt.d08','dt.dl0','dt.dl4')
    }
    out
}

k1 <- function(tt,l0,l4){
    (1-exp((l4[1]-l0[1])*tt))*l0[1]/(l4[1]-l0[1])
}
ThU_misfit <- function(tt,A,a,l0,l4,Th2U8=0,Th02i=0){
    ( 1 - exp(-l0[1]*tt)*(1-Th2U8*Th02i) - (a-1)*k1(tt,l0,l4) - A )^2
}

# converts format 1 to format 2 and vice versa
ThU_convert <- function(x){
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

get_Th230U238_ratio <- function(tt,U234U238_0){
    l4 <- lambda('U234')
    l0 <- lambda('Th230')
    a <- get_U234U238_ratio(tt,U234U238_0)
    1 - exp(-l0[1]*tt) - (a-1)*k1(tt,l0,l4)
}

get_U234U238_ratio <- function(tt,U234U238_0){
    l4 <- lambda('U234')
    1 + (U234U238_0-1)*exp(-l4[1]*tt)
}
