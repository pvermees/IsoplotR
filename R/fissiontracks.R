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
    zeta <- 4*(1+R)*MM*1e18/(Na*Lf*qap*L*dens*R)
    c(zeta,0)
}

#' Calculate the zeta calibration coefficient for fission track dating
#'
#' Determines the zeta calibration constant of a fission track dataset
#' (EDM or LA-ICP-MS) given its true age and analytical uncertainty.
#'
#' @param x an object of class \code{fissiontracks}
#' @param tst a two-element vector with the true age and its standard
#'     error
#' @param exterr logical flag indicating whether the external
#'     uncertainties associated with the age standard or the dosimeter
#'     glass (for the EDM) should be accounted for when propagating
#'     the uncertainty of the zeta calibration constant.
#' @param update logical flag indicating whether the function should
#'     return an updated version of the input data, or simply return a
#'     two-element vector with the calibration constant and its
#'     standard error.
#' @param sigdig number of significant digits
#' @return an object of class \code{fissiontracks} with an updated
#'     \code{x$zeta} value
#' @examples
#' data(examples)
#' print(examples$FT1$zeta)
#' FT <- set.zeta(examples$FT1,tst=c(250,5))
#' print(FT$zeta)
#' @export
set.zeta <- function(x,tst=c(0,0),exterr=TRUE,update=TRUE,sigdig=2){
    N <- length(x$Ns)
    L8 <- lambda('U238')[1]
    tt <- tst[1]
    st <- tst[2]
    if (!exterr) st <- 0
    if (x$format==1){
        Ns <- sum(x$x[,'Ns'])
        Ni <- sum(x$x[,'Ni'])
        rhoD <- x$rhoD
        if (!exterr) rhoD[2] <- 0
        zeta <- 2e6*(exp(L8*tt)-1)/(L8*rhoD[1]*Ns/Ni)
        zetaErr <- zeta * sqrt( (L8*exp(L8*tt)*st/(exp(L8*tt)-1))^2 +
                                (rhoD[2]/rhoD[1])^2 + 1/Ns + 1/Ni )
    } else {
        Ns <- sum(x$Ns)
        UsU <- get.UsU(x)
        UA <- sum(UsU[,1]*x$A)
        UAerr <- sqrt( sum(UsU[,2]*x$A)^2 )
        zeta <- (exp(L8*tt)-1)/(L8*Ns/(2*UA))
        zetaErr <- zeta * sqrt( ((L8*exp(L8*tt)*st)/(exp(L8*tt)-1))^2 +
                                1/Ns + (UAerr/UA)^2 )
    }
    zsz <- roundit(zeta,zetaErr,sigdig=sigdig)
    if (update){
        out <- x
        out$zeta <- c(zsz$x,zsz$err)
    } else {
        out <- matrix(c(zsz$x,zsz$err),1,2)
        colnames(out) <- c('zeta','s[zeta]')
    }
    out
}

ICP.age <- function(x,i=NA,sigdig=NA,exterr=TRUE){
    ngrains <- length(x$Ns)
    tt <- rep(NA,ngrains)
    st <- rep(NA,ngrains)
    if (exterr){
        zeta <- x$zeta
    } else {
        zeta <- c(x$zeta[1],0)
    }
    ipos <- which(x$Ns>0)
    izero <- which(x$Ns<1)
    UsU <- get.UsU(x)
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

get.UsU <- function(x){
    Aicp <- pi*(x$spotSize/2)^2
    n <- length(x$U)
    nspots <- length(unlist(x$U))
    do.average <- (nspots>n)
    out <- matrix(0,n,2)
    colnames(out) <- c('U','sU')
    m <- rep(0,n)
    if (do.average) {
        uhat <- rep(0,n)
        num <- 0
        den <- 0
    }
    for (j in 1:n){
        if (do.average){
            Uj <- stats::na.omit(x$U[[j]])
            m[j] <- length(Uj) # spots per grain
            uhat[j] <- mean(log(Uj))
            num <- num + sum((log(Uj)-uhat[j])^2)
            den <- den + m[j] - 1
        } else {
            out[j,'U'] <- x$U[[j]]
            out[j,'sU'] <- x$sU[[j]]
        }
    }
    if (do.average){
        out[,'U'] <- exp(uhat)
        vhat <- rep(num/den,n)
        for (j in 1:n){
            suhat <- x$sU[[j]]/x$U[[j]]
            vhat[j] <- vhat[j]*(1-m[j]*Aicp/x$A[j])^2 +
                       sum(suhat^2,na.rm=TRUE)*(Aicp/x$A[j])^2
        }
        out[,'sU'] <- exp(uhat)*sqrt(vhat)
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
