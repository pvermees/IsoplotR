fissiontrack.age <- function(x,i=NA,sigdig=NA,exterr=TRUE,mineral=NA){
    if (x$format < 2){
        out <- EDM.age(x,i,sigdig,exterr)
    } else if (x$format > 1){
        if (x$format == 3) {
            if (is.na(mineral)) mineral <- 'apatite'
            x$zeta <- get.absolute.zeta(mineral);
        }
        out <- ICP.age(x,i,sigdig,exterr)
    }
    out
}

get.absolute.zeta <- function(mineral){
    R <- iratio('U238U235')[1]
    MM <- imass('U')[1]
    qap <- etchfact(mineral)
    L <- tracklength(mineral)
    Lf <- lambda('fission')[1]
    dens <- mindens(mineral)
    Na <- 6.02214e23
    zeta <- MM*(1+R)*1e18/(Na*Lf*qap*L*dens*R)
    c(zeta,0)
}

ICP.age <- function(x,i=NA,sigdig=NA,exterr=TRUE){
    ngrains <- length(x$Ns)
    tt <- rep(0,ngrains)
    st <- rep(0,ngrains)
    if (exterr){
        zeta <- x$zeta
    } else {
        zeta <- c(x$zeta[1],0)
    }
    UsU <- get.UsU(x$U,x$sU)
    ipos <- which(x$Ns>0)
    izero <- which(x$Ns<1)
    # first calculate the ages of the non-zero track data:
    for (i in ipos){
        tst <- get.ICP.age(x$Ns[i],x$A[i],UsU[i,],zeta)
        tt[i] <- tst[1]
        st[i] <- tst[2]
    }
    # then use the equivalent induced track approach:
    for (i in izero){
        rho <- UsU[i,'U']/(x$A[i]*UsU[i,'sU']^2)
        Ni <- x$A[i]*UsU[i,'U']*rho
        tst <- get.EDM.age(x$Ns[i],Ni,zeta*1e6,c(rho,0))
        tt[i] <- tst[1]
        st[i] <- tst[2]
    }
    tst <- roundit(tt,st,sigdig=sigdig)
    out <- cbind(tst$x,tst$err)
    colnames(out) <- c('t','s[t]')
    out
}

get.UsU <- function(U,sU){
    ngrains <- length(U)
    nspots <- length(unlist(U))
    do.average <- (nspots>ngrains)
    out <- matrix(0,ngrains,2)
    colnames(out) <- c('U','sU')
    if (do.average) {
        muhat <- rep(0,ngrains)
        num <- 0
        den <- 0
    }
    for (i in 1:ngrains){
        if (do.average){
            sppg <- length(U[[i]]) # spots per grain
            muhat[i] <- mean(log(U[[i]]))
            num <- num + sum((log(U[[i]])-muhat[i])^2)
            den <- den + max(1,sppg-1)
        } else {
            out[i,'U'] <- U[[i]]
            out[i,'sU'] <- sU[[i]]
        }
    }
    if (do.average){
        vhat <- num/den
        out[,'U'] <- exp(muhat) #exp(muhat+vhat/2)
        out[,'sU'] <- out[,'U']*sqrt(vhat) #sqrt(exp(vhat)-1)
    }
    out
}

EDM.age <- function(x,i=NA,sigdig=2,exterr=TRUE){
    ns <- nrow(x$x)
    out <- matrix(0,ns,2)
    colnames(out) <- c('t','s[t]')
    if (exterr){
        zeta <- x$zeta
        rhoD <- x$rhoD
    } else {
        zeta <- c(x$zeta[1],0)
        rhoD <- c(x$rhoD[1],0)
    }
    for (j in 1:ns){
        tt <- get.EDM.age(x$x[j,'Ns'],x$x[j,'Ni'],zeta,rhoD)
        t.out <- roundit(tt[1],tt[2],sigdig=sigdig)
        out[j,] <- c(t.out$x,t.out$err)
    }
    if (!is.na(i)) out <- out[i,]
    out
}

# zeta and rhoD are two-element vectors
get.EDM.age <- function(Ns,Ni,zeta,rhoD){
    L8 <- lambda('U238')[1]
    if (Ns<1){
        Ns <- Ns+0.5
        Ni <- Ni+0.5
    }
    tt <- log(1+0.5*L8*(zeta[1]/1e6)*rhoD[1]*(Ns/Ni))/L8
    st <- tt*sqrt(1/Ns + 1/Ni + (rhoD[2]/rhoD[1])^2 + (zeta[2]/zeta[1])^2)
    c(tt,st)
}
# zeta is a two-element vector
get.ICP.age <- function(Ns,A,UsU,zeta){
    L8 <- lambda('U238')[1]
    tt <- log(1+L8*zeta[1]*Ns/(2*UsU[1]*A))/L8
    st <- tt * sqrt(1/Ns + (zeta[2]/zeta[1])^2 + (UsU[2]/UsU[1])^2)
    c(tt,st)
}
