# convert isotope dilution derived concentrations to ratios x = matrix
# with columns 'Luppm','errLuppm', 'Hfppm','errHfppm' and
# 'Hf176Hf177','errHf176Hf177'
ppm2ratios.LuHf <- function(x,exterr=FALSE,common=FALSE,...){
    R65Lu <- iratio('Lu176Lu175')[1]
    R47Hf <- iratio('Hf174Hf177')[1]
    R67Hf0 <- iratio('Hf176Hf177')[1]
    R87Hf <- iratio('Hf178Hf177')[1]
    R97Hf <- iratio('Hf179Hf177')[1]
    R07Hf <- iratio('Hf180Hf177')[1]
    if (common){ # make common Hf correction
        R67Hf <-  x$x[,'Hf176Hf177'] - R67Hf0
        dR67Hf.dR67Hf0 <- -1
    } else {
        R67Hf <-  x$x[,'Hf176Hf177']
        dR67Hf.dR67Hf0 <- 0
    }
    Lu <- x$x[,'Luppm']
    Hf <- x$x[,'Hfppm']
    MMLu <- imass('Lu')[1]
    MMHf <- imass('Hf')[1]
    invf6Lu <- 1 + 1/R65Lu
    invf7Hf <- R47Hf + R67Hf + 1 + R87Hf + R97Hf + R07Hf
    Lu176Hf177 <- (invf7Hf/invf6Lu)*(Lu/Hf)*(MMHf/MMLu)
    dinvf6Lu.dR65Lu <- -1/R65Lu^2
    dinvf7Hf.dR47Lu <- 1
    dinvf7Hf.dR67Hf <- 1
    dinvf7Hf.dR67Hf0 <- dinvf7Hf.dR67Hf*dR67Hf.dR67Hf0
    dinvf7Hf.dR87Hf <- 1
    dinvf7Hf.dR97Hf <- 1
    dinvf7Hf.dR07Hf <- 1
    E <- matrix(0,11,11)
    J <- matrix(0,2,11)
    if (exterr){
        E[1,1] <- iratio('Lu176Lu175')[2]^2 # var(R65Lu)
        E[2,2] <- iratio('Hf174Hf177')[2]^2 # var(R47Hf)
        E[3,3] <- iratio('Hf176Hf177')[2]^2 # var(R67Hf0)
        E[4,4] <- iratio('Hf178Hf177')[2]^2 # var(R87Hf)
        E[5,5] <- iratio('Hf179Hf177')[2]^2 # var(R97Hf)
        E[6,6] <- iratio('Hf180Hf177')[2]^2 # var(R07Hf)
        E[7,7] <- imass('Lu')[2]^2 # var(MMLu)
        E[8,8] <- imass('Hf')[2]^2 # var(MMHf)
    }
    nn <- dim(x$x)[1]
    out <- matrix(NA,nn,5) # X=Lu176/Hf177 and Y=Hf176/Hf177
    colnames(out) <- c('X','sX','Y','sY','rXY')
    out[,1] <- Lu176Hf177
    out[,3] <- R67Hf
    for (i in 1:nn){
        E[9,9] <- x$x[i,'errHf176Hf177']^2 # var(R67Hf)
        E[10,10] <- x$x[i,'errLuppm']^2    # var(Lu)
        E[11,11] <- x$x[i,'errHfppm']^2    # var(Hf)
        J[1,1] <- -Lu176Hf177[i]*dinvf6Lu.dR65Lu/invf6Lu
        J[1,2] <-  Lu176Hf177[i]*dinvf7Hf.dR47Lu/invf7Hf[i]
        J[1,3] <-  Lu176Hf177[i]*dinvf7Hf.dR67Hf0/invf7Hf[i]
        J[1,4] <-  Lu176Hf177[i]*dinvf7Hf.dR87Hf/invf7Hf[i]
        J[1,5] <-  Lu176Hf177[i]*dinvf7Hf.dR97Hf/invf7Hf[i]
        J[1,6] <-  Lu176Hf177[i]*dinvf7Hf.dR07Hf/invf7Hf[i]
        J[1,7] <- -Lu176Hf177[i]/MMLu
        J[1,8] <-  Lu176Hf177[i]/MMHf
        J[1,9] <-  Lu176Hf177[i]*dinvf7Hf.dR67Hf/invf7Hf[i]
        J[1,10] <-  Lu176Hf177[i]/Lu[i]
        J[1,11] <- -Lu176Hf177[i]/Hf[i]
        J[2,3] <- dR67Hf.dR67Hf0
        J[2,9] <- 1
        E2 <- J %*% E %*% t(J)
        if (all(is.finite(E2))){
            out[i,2] <- sqrt(E2[1,1])
            out[i,4] <- sqrt(E2[2,2])
            out[i,5] <- stats::cov2cor(E2)[1,2]
        }
    }
    out
}

get.LuHf.ratio <- function(tt,st,exterr=TRUE){
    get.PD.ratio(tt,st,nuclide='Lu176',exterr=exterr)
}

get.LuHf.age <- function(Hf176Lu176,sHf176Lu176,exterr=TRUE){
    get.PD.age(Hf176Lu176,sHf176Lu176,nuclide='Lu176',exterr=exterr)
}

LuHf.age <- function(x,exterr=TRUE,i=NA,sigdig=NA,i2i=TRUE,...){
    PD.age(x,'Lu176',exterr=exterr,i=i,sigdig=sigdig,i2i=i2i,...)
}
