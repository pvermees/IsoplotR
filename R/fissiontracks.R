fissiontrack.age <- function(x,i=NA,sigdig=NA,exterr=TRUE){
    if (x$format==1){
        out <- EDM.age(x,i,sigdig,exterr)
    } else if (x$format==2){
        out <- ICP.age(x,i,sigdig,exterr)
    }
    out
}

ICP.age <- function(x,i=NA,sigdig=NA,exterr=TRUE){
    l238 <- lambda('U238')[1]
    ngrains <- length(x$Ns)
    out <- matrix(0,ngrains,2)
    colnames(out) <- c('t','s[t]')
    if (exterr){
        zeta <- x$zeta
    } else {
        zeta <- c(x$zeta[1],0)
    }
    UsU <- get.UsU(x$U,x$sU)
    tt <- log(1+l238*(zeta[1]*100)*x$Ns/(2*UsU[,'U']*x$A))/l238
    st <- tt * sqrt(1/x$Ns + (zeta[2]/zeta[1])^2 + UsU[,'sU']/UsU[,'U'])
    tst <- roundit(tt,st,sigdig=sigdig)
    out[,'t'] <- tst$x
    out[,'s[t]'] <- tst$err
    out
}

get.UsU <- function(U,sU){
    ngrains <- length(U)
    nspots <- length(unlist(x$U))
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
        out[,'U'] <- exp(muhat+vhat/2)
        out[,'sU'] <- out[,'U']*sqrt(exp(vhat)-1)
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
    out <- c(0,0)
    out[1] <- log(1+0.5*L8*(zeta[1]/1e6)*rhoD[1]*(Ns/Ni))/L8
    out[2] <- out[1]*sqrt(1/Ns + 1/Ni +
                          (rhoD[2]/rhoD[1])^2 +
                          (zeta[2]/zeta[1])^2)
    out
}
