# convert isotope dilution derived concentrations to ratios
# x = matrix with columns 'Rbppm','errRbppm', 'Srppt','errSrppt' 
# and 'Sr87Sr86','errSr87Sr86'
ppm2ratios.RbSr <- function(x,exterr=FALSE,common=FALSE,...){
    R57Rb <- iratio('Rb85Rb87')[1]
    R46Sr <- iratio('Sr84Sr86')[1]
    R76Sr0 <- iratio('Sr87Sr86')[1]
    R86Sr <- iratio('Sr88Sr86')[1]
    if (common){ # make common Sr correction
        R76Sr <-  x$x[,'Sr87Sr86'] - R76Sr0
        dR76Sr.dR76Sr0 <- -1
    } else {
        R76Sr <- x$x[,'Sr87Sr86']
        dR76Sr.dR76Sr0 <- 0
    }
    Rb <- x$x[,'Rbppm']
    Sr <- x$x[,'Srppm']
    MMRb <- imass('Rb')[1]
    MMSr <- imass('Sr')[1]
    invf7Rb <- 1 + R57Rb
    invf6Sr <- R46Sr + 1 + R76Sr + R86Sr
    Rb87Sr86 <- (invf6Sr/invf7Rb)*(Rb/Sr)*(MMSr/MMRb)
    dinvf7Rb.dR57Rb <- 1
    dinvf6Sr.dR46Sr <- 1
    dinvf6Sr.dR76Sr <- 1
    dinvf6Sr.dR76Sr0 <- dinvf6Sr.dR76Sr*dR76Sr.dR76Sr0
    dinvf6Sr.dR86Sr <- 1
    E <- matrix(0,9,9)
    J <- matrix(0,2,9)
    if (exterr){
        E[1,1] <- iratio('Rb85Rb87')[2]^2 # var(R57Rb)
        E[2,2] <- iratio('Sr84Sr86')[2]^2 # var(R46Sr)
        E[3,3] <- iratio('Sr87Sr86')[2]^2 # var(R76Sr0)
        E[4,4] <- iratio('Sr88Sr86')[2]^2 # var(R86Sr)
        E[5,5] <- imass('Rb')[2]^2 # var(MMRb)
        E[6,6] <- imass('Sr')[2]^2 # var(MMSr)
    }
    nn <- dim(x$x)[1]
    out <- matrix(NA,nn,5) # X=Rb87/Sr86 and Y=Sr87/Sr86
    colnames(out) <- c('X','sX','Y','sY','rXY')
    out[,1] <- Rb87Sr86
    out[,3] <- R76Sr
    for (i in 1:nn){
        E[7,7] <- x$x[i,'errSr87Sr86']^2 # var(R76Sr)
        E[8,8] <- x$x[i,'errRbppm']^2
        E[9,9] <- x$x[i,'errSrppm']^2
        J[1,1] <- -Rb87Sr86[i]*dinvf7Rb.dR57Rb/invf7Rb
        J[1,2] <-  Rb87Sr86[i]*dinvf6Sr.dR46Sr/invf6Sr[i]
        J[1,3] <-  Rb87Sr86[i]*dinvf6Sr.dR76Sr0/invf6Sr[i]
        J[1,4] <-  Rb87Sr86[i]*dinvf6Sr.dR86Sr/invf6Sr[i]
        J[1,5] <- -Rb87Sr86[i]/MMRb
        J[1,6] <-  Rb87Sr86[i]/MMSr
        J[1,7] <-  Rb87Sr86[i]*dinvf6Sr.dR76Sr/invf6Sr[i]
        J[1,8] <-  Rb87Sr86[i]/Rb[i]
        J[1,9] <- -Rb87Sr86[i]/Sr[i]
        J[2,3] <-  dR76Sr.dR76Sr0
        J[2,7] <-  1
        E2 <- J %*% E %*% t(J)
        if (all(is.finite(E2))){
            out[i,2] <- sqrt(E2[1,1])
            out[i,4] <- sqrt(E2[2,2])
            out[i,5] <- stats::cov2cor(E2)[1,2]
        }
    }
    out
}

get.RbSr.ratio <- function(tt,st,exterr=TRUE){
    get.PD.ratio(tt,st,'Rb87',exterr)
}

get.RbSr.age <- function(Rb87Sr86,sRb87Sr86,exterr=TRUE){
    get.PD.age(Rb87Sr86,sRb87Sr86,'Rb87',exterr=exterr)
}

RbSr.age <- function(x,exterr=TRUE,i=NA,sigdig=NA,i2i=TRUE,...){
    PD.age(x,'Rb87',exterr=exterr,i=i,sigdig=sigdig,i2i=i2i,...)
}
