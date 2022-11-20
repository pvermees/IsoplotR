# Ludwig regression with initial disequilibrium
ludwigd <- function(x,model=1,anchor=0,bayes=FALSE){
    if (x$format<4){
        init <- init.ludwigd(x,anchor=anchor)
        fit <- optim(par=init,fn=LL.ludwigd,method="L-BFGS-B",
                     lower=c(init['tt']/10,init['a0']/10,0),
                     upper=c(init['tt']*10,init['a0']*10,20),
                     hessian=TRUE,x=x,anchor=anchor)
    } else if (x$format<7){
        # TODO
    } else {
        # TODO
    }
    if (det(fit$hessian)<0){
        warning('Ill-conditioned Hessian, replaced by ',
                'nearest positive definite matrix')
        fit$hessian <- nearPD(fit$hessian)
    }
    fit
}

init.ludwigd <- function(x,model=1,anchor=0){
    if (x$format<4){
        yd <- data2york(x,option=2)
        yfit <- york(yd)
        PbU0 <- -yfit$b[1]/yfit$a[1]
        out <- vector()
        if (anchor[1]==2 && length(anchor)>2 && anchor[3]>0){
            out['tt'] <- anchor[2]
            out['a0'] <- yfit$a[1]
        } else if (anchor[1]==1 && iratio('Pb207Pb206')[2]>0){
            out['tt'] <- get.Pb206U238.age(x=PbU0)[1]
            out['a0'] <- iratio('Pb207Pb206')[1]
        } else {
            out['tt'] <- get.Pb206U238.age(x=PbU0)[1]
            out['a0'] <- yfit$a[1]
        }
        if (model==3){
            pars['w'] <- 0.01
        }
        if (x$d$U48$option==1 && x$d$U48$sx>0){
            out['U48i'] <- x$d$U48$x
        } else if (x$d$U48$option==2){
            if (x$d$U48$sx>0){
                out['U48i'] <- 1
            } else {
                stop('Zero uncertainty of measured 234/238 activity ratio')
            }
        }
        if (x$d$ThU$option==1 && x$d$ThU$sx>0){
            out['ThUi'] <- x$d$ThU$x
        } else if (x$d$ThU$option==2){
            if (x$d$ThU$sx>0){
                out['ThUi'] <- 1
            } else {
                stop('Zero uncertainty of measured 230/238 activity ratio')
            }
        }
    }
    out
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
                LL <- LL + dnorm(x=a0,
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
                LL <- LL + dnorm(x=tt,mean=anchor[2],sd=st,log=TRUE)
            } else {
                tt <- anchor[2]
            }
        } else {
            tt <- pars['tt']
            a0 <- pars['a0']
        }
        ta0w <- c('tt'=unname(tt),'a0'=unname(a0))
        if (model==3){
            ta0w['w'] <- pars['w']
        }
        LL <- LL + data2ludwigd(X,ta0w,exterr=exterr)$LL
        if (x$d$U48$option==1 && x$d$U48$sx>0){
            U48i <- pars['U48i']
            LL <- LL + dnorm(x=U48i,mean=x$d$U48$x,sd=x$d$U48$sx,log=TRUE)
        } else if (x$d$U48$option==2 && x$d$U48$sx>0){
            pred <- mclean(tt=tt,d=X$d)
            LL <- LL + dnorm(x=pred$U48,mean=x$d$U48$x,sd=x$d$U48$sx,log=TRUE)
        }
        if (x$d$ThU$option==1 && x$d$ThU$sx>0){
            ThUi <- pars['ThUi']
            LL <- LL + dnorm(x=ThUi,mean=x$d$ThU$x,sd=x$d$ThU$sx,log=TRUE)
        } else if (x$d$ThU$option==2 && x$d$ThU$sx>0){
            pred <- mclean(tt=tt,d=X$d)
            LL <- LL + dnorm(x=pred$ThU,mean=x$d$ThU$x,sd=x$d$ThU$sx,log=TRUE)
        }
    }
    -LL
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
        b0 <- ta0b0w['w']
        Z <- zeros
        L0 <- zeros
        NP <- 3 # tt, a0, b0
        NR <- 3 # X, Y, Z
    } else if (x$format%in%c(7,8)){
        b0 <- ta0b0w['w']
        Z <- zeros
        W <- zeros
        L0 <- zeros
        NP <- 3 # tt, a0, b0
        NR <- 4 # X, Y, Z, W
    } else {
        stop('Incorrect input format.')
    }
    np <- min(NP+1,length(ta0b0w))  # np = NP (no w) or NP+1 (with w)
    if (np==(NP+1)) w <- ta0b0w[np] # model 3
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
        Ew <- get.Ew(w=w,format=x$format,ns=ns,D=D)
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
    out$LL <- -(NP*ns*log(2*pi) + detED + out$SS)/2
    out
}
