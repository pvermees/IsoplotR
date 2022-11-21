# Ludwig regression with initial disequilibrium
ludwigd <- function(x,model=1,anchor=0,exterr=FALSE,bayes=FALSE){
    init <- init.ludwigd(x,model=model,anchor=anchor)
    fit <- optim(par=init$pars,fn=LL.ludwigd,method="L-BFGS-B",
                 lower=init$lower,upper=init$upper,hessian=TRUE,
                 x=x,anchor=anchor,model=model,exterr=exterr)
    if (det(fit$hessian)<0){
        warning('Ill-conditioned Hessian, replaced by ',
                'nearest positive definite matrix')
        fit$hessian <- nearPD(fit$hessian)
    }
    fit
}

init.ludwigd <- function(x,model=1,anchor=0){
    pars <- lower <- upper <- vector()
    if (x$format<4){
        yd <- data2york(x,option=2)
        yfit <- york(yd)
        PbU0 <- -yfit$b[1]/yfit$a[1]
        if (anchor[1]==1){
            pars['tt'] <- get.Pb206U238.age(x=PbU0)[1]
            lower['tt'] <- pars['tt']/10
            upper['tt'] <- pars['tt']*10
            if (iratio('Pb207Pb206')[2]>0){
                a0 <- iratio('Pb207Pb206')[1]
                sa0 <- iratio('Pb207Pb206')[2]
                pars['a0'] <- a0
                lower['a0'] <- max(0,a0-sa0*3)
                upper['a0'] <- a0+sa0*3
            }
        } else if (anchor[1]==2){
            if (length(anchor)>2 && anchor[3]>0){
                pars['tt'] <- anchor[2]
                lower['tt'] <- pars['tt']-anchor[3]*3
                upper['tt'] <- pars['tt']+anchor[3]*3
            }
            pars['a0'] <- yfit$a[1]
            lower['a0'] <- pars['a0']/10
            upper['a0'] <- pars['a0']*10
        } else {
            pars['tt'] <- get.Pb206U238.age(x=PbU0)[1]
            lower['tt'] <- pars['tt']/10
            upper['tt'] <- pars['tt']*10            
            pars['a0'] <- yfit$a[1]
            lower['a0'] <- pars['a0']/10
            upper['a0'] <- pars['a0']*10            
        }
    } else if (x$format<7){
        yda <- data2york(x,option=3)
        yfita <- york(yda)
        Pb6U8 <- -yfita$b[1]/yfita$a[1]
        ydb <- data2york(x,option=4)
        yfitb <- york(ydb)
        if (anchor[1]==1){
            pars['tt'] <- get.Pb206U238.age(x=Pb6U8)[1]
            lower['tt'] <- pars['tt']/10
            upper['tt'] <- pars['tt']*10
            if (iratio('Pb206Pb204')[2]>0){
                a0 <- iratio('Pb206Pb204')[1]
                sa0 <- iratio('Pb206Pb204')[2]
                pars['a0'] <- a0
                lower['a0'] <- max(0,a0-sa0*3)
                upper['a0'] <- a0+sa0*3
            }
            if (iratio('Pb207Pb204')[2]>0){
                b0 <- iratio('Pb207Pb204')[1]
                sb0 <- iratio('Pb207Pb204')[2]
                pars['b0'] <- b0
                lower['b0'] <- max(0,b0-sb0*3)
                upper['b0'] <- b0+sb0*3
            }
        } else if (anchor[1]==2){
            if (length(anchor)>2 && anchor[3]>0){
                pars['tt'] <- anchor[2]
                lower['tt'] <- pars['tt']-anchor[3]*3
                upper['tt'] <- pars['tt']+anchor[3]*3
            }
            pars['a0'] <- 1/yfita$a[1]
            lower['a0'] <- pars['a0']/10
            upper['a0'] <- pars['a0']*10
            pars['b0'] <- 1/yfitb$a[1]
            lower['b0'] <- pars['b0']/10
            upper['b0'] <- pars['b0']*10
        } else {
            pars['tt'] <- get.Pb206U238.age(x=Pb6U8)[1]
            lower['tt'] <- pars['tt']/10
            upper['tt'] <- pars['tt']*10
            pars['a0'] <- 1/yfita$a[1]
            lower['a0'] <- pars['a0']/10
            upper['a0'] <- pars['a0']*10
            pars['b0'] <- 1/yfitb$a[1]
            lower['b0'] <- pars['b0']/10
            upper['b0'] <- pars['b0']*10
        }
    } else {
        yd <- data2york(x,option=2)
        yfit <- york(yd)
        Pb6U8 <- -yfit$b[1]/yfit$a[1]
        tt <- get.Pb206U238.age(x=Pb6U8)[1]
        yda <- data2york(x,option=6,tt=tt)
        yfita <- york(yda)
        ydb <- data2york(x,option=7,tt=tt)
        yfitb <- york(ydb)
        if (anchor[1]==1){
            pars['tt'] <- tt
            lower['tt'] <- tt/10
            upper['tt'] <- tt*10
            if (iratio('Pb208Pb206')[2]>0){
                a0 <- 1/iratio('Pb208Pb206')[1]
                sa0 <- iratio('Pb208Pb206')[2]*a0^2
                pars['a0'] <- a0
                lower['a0'] <- max(0,a0-sa0*3)
                upper['a0'] <- a0+sa0*3
            }
            if (iratio('Pb207Pb206')[2]>0){
                b0 <- 1/iratio('Pb208Pb207')[1]
                sb0 <- iratio('Pb208Pb207')[2]/b0^2
                pars['b0'] <- b0
                lower['b0'] <- max(0,b0-sb0*3)
                upper['b0'] <- b0+sb0*3
            }
        } else if (anchor[1]==2){
            if (length(anchor)>2 && anchor[3]>0){
                pars['tt'] <- anchor[2]
                lower['tt'] <- pars['tt']-anchor[3]*3
                upper['tt'] <- pars['tt']+anchor[3]*3
            }
            pars['a0'] <- 1/yfita$a[1]
            lower['a0'] <- pars['a0']/10
            upper['a0'] <- pars['a0']*10
            pars['b0'] <- 1/yfitb$a[1]
            lower['b0'] <- pars['b0']/10
            upper['b0'] <- pars['b0']*10
        } else {
            pars['tt'] <- tt
            lower['tt'] <- tt/10
            upper['tt'] <- tt*10
            pars['a0'] <- 1/yfita$a[1]
            lower['a0'] <- pars['a0']/10
            upper['a0'] <- pars['a0']*10
            pars['b0'] <- 1/yfitb$a[1]
            lower['b0'] <- pars['b0']/10
            upper['b0'] <- pars['b0']*10
        }
    }
    if (model==3){
        pars['w'] <- 0.01
        lower['w'] <- 0
        upper['w'] <- 10
    }
    if (x$d$U48$option==1 && x$d$U48$sx>0){
        pars['U48i'] <- x$d$U48$x
        lower['U48i'] <- max(0,pars['U48i']-x$d$U48$sx*3)
        upper['U48i'] <- pars['U48i']+x$d$U48$sx*3
    } else if (x$d$U48$option==2){
        if (x$d$U48$sx>0){
            pars['U48i'] <- 1
            lower['U48i'] <- 0
            upper['U48i'] <- 20
        } else {
            stop('Zero uncertainty of measured 234/238 activity ratio')
        }
    }
    if (x$d$ThU$option==1 && x$d$ThU$sx>0){
        pars['ThUi'] <- x$d$ThU$x
        lower['ThUi'] <- max(0,pars['ThUi']-x$d$ThU$sx*3)
        upper['ThUi'] <- pars['ThUi']+x$d$ThU$sx*3
    } else if (x$d$ThU$option==2){
        if (x$d$ThU$sx>0){
            pars['ThUi'] <- 1
            lower['ThUi'] <- 0
            upper['ThUi'] <- 20
        } else {
            stop('Zero uncertainty of measured 230/238 activity ratio')
        }
    }
    list(pars=pars,lower=lower,upper=upper)
}

LL.ludwigd <- function(pars,x,model=1,exterr=FALSE,anchor=0){
    X <- x
    if ('U48i' %in% names(pars)){
        X$d$U48$x <- pars['U48i']
        X$d$U48$option <- 1
    }
    if ('logit[ThUi]' %in% names(pars)){
        X$d$ThU$x <- pars['ThUi']
        X$d$ThU$option <- 1
    }
    LL <- 0
    if (x$format<4){
        if (anchor[1]==1){
            tt <- pars['tt']
            if (iratio('Pb207Pb206')[2]>0){
                a0 <- pars['a0']
                LL <- LL - dnorm(x=a0,
                                 mean=iratio('Pb207Pb206')[1],
                                 sd=iratio('Pb207Pb206')[2],
                                 log=TRUE)
            } else {
                a0 <- iratio('Pb207Pb206')[1]
            }
        } else if (anchor[1]==2){
            a0 <- pars['a0']
            if (length(anchor)>2 && anchor[3]>0){
                tt <- pars['tt']
                st <- anchor[3]
                LL <- LL - dnorm(x=tt,mean=anchor[2],sd=st,log=TRUE)
            } else {
                tt <- anchor[2]
            }
        } else {
            tt <- pars['tt']
            a0 <- pars['a0']
        }
        ta0b0w <- c('tt'=unname(tt),'a0'=unname(a0))
    } else if (x$format<7){
        if (anchor[1]==1){
            tt <- pars['tt']
            if (iratio('Pb206Pb204')[2]>0){
                a0 <- pars['a0']
                LL <- LL - dnorm(x=a0,
                                 mean=iratio('Pb206Pb204')[1],
                                 sd=iratio('Pb206Pb204')[2],
                                 log=TRUE)
            } else {
                a0 <- iratio('Pb206Pb204')[1]
            }
            if (iratio('Pb207Pb204')[2]>0){
                b0 <- pars['b0']
                LL <- LL - dnorm(x=b0,
                                 mean=iratio('Pb207Pb204')[1],
                                 sd=iratio('Pb207Pb204')[2],
                                 log=TRUE)
            } else {
                b0 <- iratio('Pb207Pb204')[1]
            }
        } else if (anchor[1]==2){
            a0 <- pars['a0']
            b0 <- pars['b0']
            if (length(anchor)>2 && anchor[3]>0){
                tt <- pars['tt']
                st <- anchor[3]
                LL <- LL - dnorm(x=tt,mean=anchor[2],sd=st,log=TRUE)
            } else {
                tt <- anchor[2]
            }
        } else {
            tt <- pars['tt']
            a0 <- pars['a0']
            b0 <- pars['b0']
        }
        ta0b0w <- c('tt'=unname(tt),'a0'=unname(a0),'b0'=unname(b0))
    } else {
        if (anchor[1]==1){
            tt <- pars['tt']
            if (iratio('Pb208Pb206')[2]>0){
                a0 <- pars['a0']
                LL <- LL - dnorm(x=1/a0,
                                 mean=iratio('Pb208Pb206')[1],
                                 sd=iratio('Pb208Pb206')[2],
                                 log=TRUE)
            } else {
                a0 <- 1/iratio('Pb208Pb206')[1]
            }
            if (iratio('Pb208Pb207')[2]>0){
                b0 <- pars['b0']
                LL <- LL - dnorm(x=1/b0,
                                 mean=iratio('Pb208Pb207')[1],
                                 sd=iratio('Pb208Pb207')[2],
                                 log=TRUE)
            } else {
                b0 <- 1/iratio('Pb208Pb207')[1]
            }
        } else if (anchor[1]==2){
            a0 <- pars['a0']
            b0 <- pars['b0']
            if (length(anchor)>2 && anchor[3]>0){
                tt <- pars['tt']
                st <- anchor[3]
                LL <- LL - dnorm(x=tt,mean=anchor[2],sd=st,log=TRUE)
            } else {
                tt <- anchor[2]
            }
        } else {
            tt <- pars['tt']
            a0 <- pars['a0']
            b0 <- pars['b0']
        }
        ta0b0w <- c('tt'=unname(tt),'a0'=unname(a0),'b0'=unname(b0))
    }
    if (model==3){
        ta0b0w['w'] <- pars['w']
    }
    if (model==2){
        LL <- LL + LL.ludwig.model2(ta0b0w,X)
    } else {
        LL <- LL + data2ludwigd(X,ta0b0w,exterr=exterr)$LL
    }
    if (x$d$U48$option==1 && x$d$U48$sx>0){
        U48i <- pars['U48i']
        LL <- LL - dnorm(x=U48i,mean=x$d$U48$x,sd=x$d$U48$sx,log=TRUE)
    } else if (x$d$U48$option==2 && x$d$U48$sx>0){
        pred <- mclean(tt=tt,d=X$d)
        LL <- LL - dnorm(x=pred$U48,mean=x$d$U48$x,sd=x$d$U48$sx,log=TRUE)
    }
    if (x$d$ThU$option==1 && x$d$ThU$sx>0){
        ThUi <- pars['ThUi']
        LL <- LL - dnorm(x=ThUi,mean=x$d$ThU$x,sd=x$d$ThU$sx,log=TRUE)
    } else if (x$d$ThU$option==2 && x$d$ThU$sx>0){
        pred <- mclean(tt=tt,d=X$d)
        LL <- LL - dnorm(x=pred$ThU,mean=x$d$ThU$x,sd=x$d$ThU$sx,log=TRUE)
    }
    LL
}

data2ludwigd <- function(x,ta0b0w,exterr=FALSE){
    out <- list()
    U <- iratio('U238U235')[1]
    tt <- ta0b0w['tt']
    a0 <- ta0b0w['a0']
    ns <- length(x)
    zeros <- rep(0,ns)
    X <- zeros
    Y <- zeros
    K0 <- zeros
    D <- mclean(tt=tt,d=x$d,exterr=exterr)
    if (x$format%in%c(1,2,3)){
        NP <- 2 # number of fit parameters (tt, a0)
        NR <- 2 # number of isotopic ratios (X, Y)
    } else if (x$format%in%c(4,5,6)){
        b0 <- ta0b0w['b0']
        Z <- zeros
        L0 <- zeros
        NP <- 3 # tt, a0, b0
        NR <- 3 # X, Y, Z
    } else if (x$format%in%c(7,8)){
        b0 <- ta0b0w['b0']
        Z <- zeros
        W <- zeros
        L0 <- zeros
        NP <- 3 # tt, a0, b0
        NR <- 4 # X, Y, Z, W
    } else {
        stop('Incorrect input format.')
    }
    np <- min(NP+1,length(ta0b0w))  # np = NP (no w) or NP+1 (with w)
    if (np==(NP+1)) w <- ta0b0w['w'] # model 3
    E <- matrix(0,NR*ns+7,NR*ns+7)
    J <- matrix(0,NP*ns,NR*ns+7)
    J[1:(NP*ns),1:(NP*ns)] <- diag(NP*ns)
    nc <- length(D$ThUi) # nc>1 if each aliquot has its own diseq correction
    j <- 1
    for (i in 1:ns){
        wd <- wetherill(x,i=i)
        X[i] <- wd$x['Pb207U235']
        Y[i] <- wd$x['Pb206U238']
        if (x$format%in%c(4,5,6)){
            Z[i] <- wd$x['Pb204U238']
        } else if (x$format>6){
            Z[i] <- wd$x['Pb208Th232']
            W[i] <- wd$x['Th232U238']
        }
        E[(0:(NR-1))*ns+i,(0:(NR-1))*ns+i] <- wd$cov
        if (nc>1) j <- i
        J[i,NR*ns+2] <- -D$dPb207U235dl35[j]     #dKdl35
        J[i,NR*ns+5] <- -D$dPb207U235dl31[j]     #dKdl31
        J[ns+i,NR*ns+1] <- -D$dPb206U238dl38[j]  #dLdl38
        J[ns+i,NR*ns+3] <- -D$dPb206U238dl34[j]  #dLdl34
        J[ns+i,NR*ns+6] <- -D$dPb206U238dl30[j]  #dLdl30
        J[ns+i,NR*ns+7] <- -D$dPb206U238dl26[j]  #dLdl26
        if (x$format>6) J[2*ns+i,NR*ns+4] <- -D$dPb208Th232dl32[j] # dMdl32
    }
    E[NR*ns+1:7,NR*ns+1:7] <- getEl()
    ED <- J%*%E%*%t(J)
    if (np==(NP+1)){ # fit overdispersion
        Ew <- get.Ewd(w=w,format=x$format,ns=ns,D=D)
        ED <- ED + Ew
    }
    i1 <- 1:ns
    i2 <- (ns+1):(2*ns)
    if (x$format<4){
        O <- blockinverse(AA=ED[i1,i1],BB=ED[i1,i2],
                          CC=ED[i2,i1],DD=ED[i2,i2],doall=TRUE)
    } else {
        i3 <- (2*ns+1):(3*ns)
        O <- blockinverse3x3(AA=ED[i1,i1],BB=ED[i1,i2],CC=ED[i1,i3],
                             DD=ED[i2,i1],EE=ED[i2,i2],FF=ED[i2,i3],
                             GG=ED[i3,i1],HH=ED[i3,i2],II=ED[i3,i3])
    }
    if (x$format%in%c(1,2,3)){
        K0 <- X - D$Pb207U235 + a0*U*(D$Pb206U238 - Y)
        A <- t(K0%*%(O[i1,i1]+t(O[i1,i1]))*a0*U +
               K0%*%(O[i1,i2]+t(O[i2,i1])))
        B <- -(a0*U*(O[i1,i1]+t(O[i1,i1]))*a0*U +
               (O[i2,i2]+t(O[i2,i2])) +
               a0*U*(O[i1,i2]+t(O[i1,i2])) +
               (O[i2,i1]+t(O[i2,i1]))*a0*U)
        L <- as.vector(solve(B,A))
        c0 <- Y - D$Pb206U238 - L
        K <- X - D$Pb207U235 - a0*U*c0
        KLM <- c(K,L)
    } else if (x$format%in%c(4,5,6)){
        K0 <- X - D$Pb207U235 - U*b0*Z
        L0 <- Y - D$Pb206U238 - a0*Z
        V <- t(K0%*%(O[i1,i1]+t(O[i1,i1]))*U*b0 +
               L0%*%(O[i1,i2]+t(O[i2,i1]))*U*b0 +
               K0%*%(O[i1,i2]+t(O[i2,i1]))*a0 +
               L0%*%(O[i2,i2]+t(O[i2,i2]))*a0 +
               K0%*%(O[i1,i3]+t(O[i3,i1])) +
               L0%*%(O[i2,i3]+t(O[i3,i2])))
        W <- -(U*b0*(O[i1,i1]+t(O[i1,i1]))*U*b0 +
               U*b0*(O[i1,i2]+t(O[i1,i2]))*a0 +
               U*b0*(O[i1,i3]+t(O[i1,i3])) +
               a0*(O[i2,i1]+t(O[i2,i1]))*U*b0 +
               a0*(O[i2,i2]+t(O[i2,i2]))*a0 +
               a0*(O[i2,i3]+t(O[i2,i3])) +
               (O[i3,i1]+t(O[i3,i1]))*U*b0 +
               (O[i3,i2]+t(O[i3,i2]))*a0 +
               (O[i3,i3]+t(O[i3,i3])))
        M <- as.vector(solve(W,V))
        c0 <- as.vector(Z - M)
        K <- as.vector(X - D$Pb207U235 - U*b0*c0)
        L <- as.vector(Y - D$Pb206U238 - a0*c0)
        KLM <- c(K,L,M)
    } else if (x$format%in%c(7,8)){
        Wd <- diag(W)
        K0 <- X - D$Pb207U235 - (Z-D$Pb208Th232)*U*W*b0
        L0 <- Y - D$Pb206U238 - (Z-D$Pb208Th232)*W*a0
        AA <- (Wd%*%O[i1,i1]%*%Wd)*(U*b0)^2 +
            (Wd%*%O[i2,i2]%*%Wd)*a0^2 + O[i3,i3] +
            U*a0*b0*Wd%*%(O[i1,i2]+O[i2,i1])%*%Wd +
            U*b0*(Wd%*%O[i1,i3]+O[i3,i1]%*%Wd) +
            a0*(Wd%*%O[i2,i3]+O[i3,i2]%*%Wd)
        BT <- t(U*b0*K0%*%O[i1,i1]%*%Wd +
                a0*L0%*%O[i2,i2]%*%Wd +
                a0*K0%*%O[i1,i2]%*%Wd +
                U*b0*L0%*%O[i2,i1]%*%Wd +
                K0%*%O[i1,i3] + L0%*%O[i2,i3])
        CC <- U*b0*Wd%*%O[i1,i1]%*%K0 +
            a0*Wd%*%O[i2,i2]%*%L0 +
            a0*Wd%*%O[i2,i1]%*%K0 +
            U*b0*Wd%*%O[i1,i2]%*%L0 +
            O[i3,i1]%*%K0 + O[i3,i2]%*%L0
        M <- as.vector(solve(-(AA+t(AA)),(BT+CC)))
        c0 <- as.vector(Z - D$Pb208Th232 - M)
        K <- as.vector(X - D$Pb207U235 - c0*b0*U*W)
        L <- as.vector(Y - D$Pb206U238 - c0*a0*W)
        KLM <- c(K,L,M)
    }
    out$c0 <- c0
    out$SS <- KLM%*%O%*%KLM
    detED <- determinant(ED,logarithm=TRUE)$modulus
    out$LL <- (NP*ns*log(2*pi) + detED + out$SS)/2
    out
}

get.Ewd <- function(w=0,format=1,ns=1,D=mclean()){
    i1 <- 1:ns
    i2 <- (ns+1):(2*ns)
    if (format<4) {
        J <- matrix(0,2*ns,ns)
    } else {
        J <- matrix(0,3*ns,ns)
        i3 <- (2*ns+1):(3*ns)
    }
    diag(J[i1,i1]) <- -D$dPb207U235dt      # dKdt
    diag(J[i2,i1]) <- -D$dPb206U238dt      # dLdt
    if (format>6) {
        diag(J[i3,i1]) <- -D$dPb208Th232dt # dMdt
    }
    dEdx <- w^2
    dEdx*J%*%t(J)
}

LL.ludwig.model2 <- function(ta0b0,x,exterr=FALSE){
    tt <- ta0b0['tt']
    a0 <- ta0b0['a0']
    nn <- length(x)
    if (x$format<4){
        xy <- data2york(x,option=2)[,c('X','Y'),drop=FALSE]
        xr <- age_to_U238Pb206_ratio(tt,st=0,d=x$d)[1]
        yr <- age_to_Pb207Pb206_ratio(tt,st=0,d=x$d)[1]
        a <- a0
        b <- (yr-a0)/xr
        # sum of the squared Deming distances:
        SS <- sum(((xy[,'Y']-b*xy[,'X']-a)^2)/(1+b^2))
        if (exterr){
            D <- mclean(tt,d=x$d,exterr=exterr)
            dypdxr <- (a0-yr)*xy[,'X',drop=FALSE]/xr^2
            dypdyr <- xy[,'X',drop=FALSE]/xr
            dxrdPbU <- -xr^2
            dyrdPbPb <- 1
            dPbUdl <- rep(0,7)
            dPbPbdl <- rep(0,7)
            dPbUdl[1] <- D$dPb206U238dl38
            dPbUdl[3] <- D$dPb206U238dl34
            dPbUdl[6] <- D$dPb206U238dl30
            dPbUdl[7] <- D$dPb206U238dl26
            dPbPbdl[1] <- D$dPb207Pb206dl38
            dPbPbdl[2] <- D$dPb207Pb206dl35
            dPbPbdl[3] <- D$dPb207Pb206dl34
            dPbPbdl[5] <- D$dPb207Pb206dl31
            dPbPbdl[6] <- D$dPb207Pb206dl30
            dPbPbdl[7] <- D$dPb207Pb206dl26
            J <- (dypdxr*dxrdPbU)%*%dPbUdl + (dypdyr*dyrdPbPb)%*%dPbPbdl
            covmat <- diag(SS/(nn-2),nn,nn) + J%*%getEl()%*%t(J)
            LL <- LL.norm(dy,covmat)
        } else {
            LL <- SS2LL(SS=SS,nn=nn)
        }
    } else {
        b0 <- ta0b0['b0']
        if (x$format<7) xy <- get.UPb.isochron.ratios.204(x)
        else xy <- get.UPb.isochron.ratios.208(x,tt=tt)[,1:4]
        x6 <- xy[,1,drop=FALSE] # U238Pb206
        y6 <- xy[,2,drop=FALSE] # Pb204Pb206 or Pb208cPb206
        x7 <- xy[,3,drop=FALSE] # U235Pb207
        y7 <- xy[,4,drop=FALSE] # Pb204Pb207 or Pb208cPb207
        r86 <- age_to_U238Pb206_ratio(tt,st=0,d=x$d)[1]
        r57 <- age_to_U235Pb207_ratio(tt,st=0,d=x$d)[1]
        y6p <- (r86-x6)/(a0*r86)
        y7p <- (r57-x7)/(b0*r57)
        dy <- cbind(y6-y6p,y7-y7p)
        E <- stats::cov(dy)
        if (exterr){
            D <- mclean(tt=tt,d=x$d,exterr=exterr)
            dy6pd68 <- -x6/a0
            dy7pd75 <- -x7/b0
            d68dl <- matrix(0,1,7)
            d68dl[1] <- D$dPb206U238dl38
            d68dl[3] <- D$dPb206U238dl34
            d68dl[6] <- D$dPb206U238dl30
            d68dl[7] <- D$dPb206U238dl26
            d75dl <- matrix(0,1,7)
            d75dl[2] <- D$dPb207U235dl35
            d75dl[5] <- D$dPb207U235dl31
            J <- matrix(0,2*nn,7)
            J[1:nn,] <- dy6pd68%*%d68dl
            J[nn+(1:nn),] <- dy7pd75%*%d75dl
            EE <- matrix(0,2*nn,2*nn)
            diag(EE)[1:nn] <- E[1,1]
            diag(EE)[nn+(1:nn)] <- E[2,2]
            diag(EE[1:nn,nn+(1:nn)]) <- E[1,2]
            diag(EE[nn+(1:nn),1:nn]) <- E[2,1]
            covmat <- EE + J %*% getEl() %*% t(J)
            LL <- LL.norm(c(y6-y6p,y7-y7p),covmat)
        } else {
            LL <- sum(apply(dy,1,LL.norm,covmat=E))
        }
    }
    LL
}
