# convert isotope dilution derived concentrations to ratios
# x = matrix with columns 'Sm[ppm]','errSm[ppm]', 'Nd[ppm]','errNd[ppm]' 
# and 'Nd143Nd144','errNd143Nd144'
ppm2ratios.SmNd <- function(x,exterr=FALSE,common=FALSE,...){
    R4452Sm <- iratio('Sm144Sm152')[1]
    R4752Sm <- iratio('Sm147Sm152')[1]
    R4852Sm <- iratio('Sm148Sm152')[1]
    R4952Sm <- iratio('Sm149Sm152')[1]
    R5052Sm <- iratio('Sm150Sm152')[1]
    R5452Sm <- iratio('Sm154Sm152')[1]
    R24Nd <- iratio('Nd142Nd144')[1]
    R34Nd0 <- iratio('Nd143Nd144')[1]
    R54Nd <- iratio('Nd145Nd144')[1]
    R64Nd <- iratio('Nd146Nd144')[1]
    R84Nd <- iratio('Nd148Nd144')[1]
    R04Nd <- iratio('Nd150Nd144')[1]
    if (common){ # make common Nd correction
        R34Nd <-  x$x[,'Nd143Nd144'] - R34Nd0
        dR34Nd.dR34Nd0 <- -1
    } else {
        R34Nd <- x$x[,'Nd143Nd144']
        dR34Nd.dR34Nd0 <- 0
    }
    Sm <- x$x[,'Smppm']
    Nd <- x$x[,'Ndppm']
    MMSm <- imass('Sm')[1]
    MMNd <- imass('Nd')[1]
    invf147Sm <- R4452Sm/R4752Sm + 1 + R4852Sm/R4752Sm +
        R4952Sm/R4752Sm + R5052Sm/R4752Sm + 1/R4752Sm + R5452Sm/R4752Sm
    invf144Nd <- R24Nd + R34Nd + 1 + R54Nd + R64Nd + R84Nd + R04Nd
    Sm147Nd144 <- (invf144Nd/invf147Sm)*(Sm/Nd)*(MMNd/MMSm)
    dinvf147Sm.dR4452Sm <- 1/R4752Sm
    dinvf147Sm.dR4752Sm <- -(R4452Sm + R4852Sm + R4952Sm + R5052Sm +
                             1 + R5452Sm)/R4752Sm^2
    dinvf147Sm.dR4852Sm <- 1/R4752Sm
    dinvf147Sm.dR4952Sm <- 1/R4752Sm
    dinvf147Sm.dR5052Sm <- 1/R4752Sm
    dinvf147Sm.dR5452Sm <- 1/R4752Sm
    dinvf144Nd.dR24Nd <- 1
    dinvf144Nd.dR34Nd <- 1
    dinvf144Nd.dR34Nd0 <- dinvf144Nd.dR34Nd*dR34Nd.dR34Nd0
    dinvf144Nd.dR54Nd <- 1
    dinvf144Nd.dR64Nd <- 1
    dinvf144Nd.dR84Nd <- 1
    dinvf144Nd.dR04Nd <- 1
    E <- matrix(0,17,17)
    J <- matrix(0,2,17)
    if (exterr){
        E[1,1] <- iratio('Sm144Sm152')[2]^2 # var(R4452Sm)
        E[2,2] <- iratio('Sm147Sm152')[2]^2 # var(R4752Sm)
        E[3,3] <- iratio('Sm148Sm152')[2]^2 # var(R4852Sm)
        E[4,4] <- iratio('Sm149Sm152')[2]^2 # var(R4952Sm)
        E[5,5] <- iratio('Sm150Sm152')[2]^2 # var(R5052Sm)
        E[6,6] <- iratio('Sm154Sm152')[2]^2 # var(R5452Sm)
        E[7,7] <- iratio('Nd142Nd144')[2]^2 # var(R24Nd)
        E[8,8] <- iratio('Nd143Nd144')[2]^2 # var(R34Nd0)
        E[9,9] <- iratio('Nd145Nd144')[2]^2 # var(R54Nd)
        E[10,10] <- iratio('Nd146Nd144')[2]^2 # var(R64Nd)
        E[11,11] <- iratio('Nd148Nd144')[2]^2 # var(R84Nd)
        E[12,12] <- iratio('Nd150Nd144')[2]^2 # var(R04Nd)
        E[13,13] <- imass('Sm')[2]^2 # var(MMSm)
        E[14,14] <- imass('Nd')[2]^2 # var(MMNd)
    }
    nn <- dim(x$x)[1]
    out <- matrix(NA,nn,5) # X=Sm147/Nd144 and Y=Nd143/Nd144
    colnames(out) <- c('X','sX','Y','sY','rXY')
    out[,1] <- Sm147Nd144
    out[,3] <- R34Nd
    for (i in 1:nn){
        E[15,15] <- x$x[i,'errNd143Nd144']^2 # var(R34Nd)
        E[16,16] <- x$x[i,'errSmppm']^2 # var(Sm)
        E[17,17] <- x$x[i,'errNdppm']^2 # var(Nd)
        J[1,1] <- -Sm147Nd144[i]*dinvf147Sm.dR4452Sm/invf147Sm
        J[1,2] <- -Sm147Nd144[i]*dinvf147Sm.dR4752Sm/invf147Sm
        J[1,3] <- -Sm147Nd144[i]*dinvf147Sm.dR4852Sm/invf147Sm
        J[1,4] <- -Sm147Nd144[i]*dinvf147Sm.dR4952Sm/invf147Sm
        J[1,5] <- -Sm147Nd144[i]*dinvf147Sm.dR5052Sm/invf147Sm
        J[1,6] <- -Sm147Nd144[i]*dinvf147Sm.dR5452Sm/invf147Sm
        J[1,7] <-  Sm147Nd144[i]*dinvf144Nd.dR24Nd/invf144Nd[i]
        J[1,8] <-  Sm147Nd144[i]*dinvf144Nd.dR34Nd0/invf144Nd[i]
        J[1,9] <-  Sm147Nd144[i]*dinvf144Nd.dR54Nd/invf144Nd[i]
        J[1,10] <- Sm147Nd144[i]*dinvf144Nd.dR64Nd/invf144Nd[i]
        J[1,11] <- Sm147Nd144[i]*dinvf144Nd.dR84Nd/invf144Nd[i]
        J[1,12] <- Sm147Nd144[i]*dinvf144Nd.dR04Nd/invf144Nd[i]
        J[1,13] <- -Sm147Nd144[i]/MMSm
        J[1,14] <-  Sm147Nd144[i]/MMNd
        J[1,15] <-  Sm147Nd144[i]*dinvf144Nd.dR34Nd/invf144Nd[i]
        J[1,16] <-  Sm147Nd144[i]/Sm[i]
        J[1,17] <- -Sm147Nd144[i]/Nd[i]
        J[2,8] <- dR34Nd.dR34Nd0
        J[2,15] <- 1
        E2 <- J %*% E %*% t(J)
        if (all(is.finite(E2))){
            out[i,2] <- sqrt(E2[1,1])
            out[i,4] <- sqrt(E2[2,2])
            out[i,5] <- stats::cov2cor(E2)[1,2]
        }
    }
    out
}

get.SmNd.ratio <- function(tt,st,exterr=TRUE){
    get.PD.ratio(tt,st,'Sm147',exterr=exterr)
}

get.SmNd.age <- function(Nd143Sm147,sNd143Sm147,exterr=TRUE){
    get.PD.age(Nd143Sm147,sNd143Sm147,'Sm147',exterr=exterr)
}

SmNd.age <- function(x,exterr=TRUE,i=NA,sigdig=NA,i2i=TRUE,...){
    PD.age(x,'Sm147',exterr=exterr,i=i,sigdig=sigdig,i2i=i2i,...)
}
