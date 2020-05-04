# convert isotope dilution derived concentrations to ratios
# x = matrix with columns 'Reppm','errReppm', 'Osppt','errOsppt' 
# and 'Os187Os188','errOs187Os188'
ppm2ratios.ReOs <- function(x,exterr=FALSE,common=FALSE,...){
    R57Re <- iratio('Re185Re187')[1]
    R42Os <- iratio('Os184Os192')[1]
    R62Os <- iratio('Os186Os192')[1]
    R72Os <- iratio('Os187Os192')[1]
    R82Os <- iratio('Os188Os192')[1]
    R92Os <- iratio('Os189Os192')[1]
    R02Os <- iratio('Os190Os192')[1]
    dR78Os.dR78Os <- 1
    if (common) { # make common Os correction
        R78Os <-  x$x[,'Os187Os188'] - R72Os/R82Os
        dR78Os.dR72Os <- -1/R82Os
        dR78Os.dR82Os <- R72Os/R82Os^2 
    } else {
        R78Os <- x$x[,'Os187Os188']
        dR78Os.dR72Os <- 0
        dR78Os.dR82Os <- 0
    }
    Re <- x$x[,'Reppm']
    Os <- x$x[,'Osppm']
    MMOs <- imass('Os')[1]
    MMRe <- imass('Re')[1]
    invf187Re <- 1+R57Re
    invf188Os <- R42Os/R82Os + R62Os/R82Os + R78Os + 1 +
                 R92Os/R82Os + R02Os/R82Os + 1/R82Os
    Re187Os188 <- (invf188Os/invf187Re)*(Re/Os)*(MMOs/MMRe)
    dinvf187Re.dR57Re <- 1
    dinvf188Os.dR42Os <- 1/R82Os
    dinvf188Os.dR62Os <- 1/R82Os
    dinvf188Os.dR72Os <- 1
    dinvf188Os.dR82Os <- -(R42Os + R82Os + R92Os + 1)/R82Os^2
    dinvf188Os.dR92Os <- 1/R82Os
    dinvf188Os.dR02Os <- 1/R82Os
    dinvf188Os.dR78Os <- 1
    E <- matrix(0,12,12)
    J <- matrix(0,2,12)
    if (exterr){
        E[1,1] <- iratio('Re185Re187')[2]^2 # var(R57Re)
        E[2,2] <- iratio('Os184Os192')[2]^2 # var(R42Os)
        E[3,3] <- iratio('Os186Os192')[2]^2 # var(R62Os)
        E[4,4] <- iratio('Os187Os192')[2]^2 # var(R72Os)
        E[5,5] <- iratio('Os188Os192')[2]^2 # var(R82Os)
        E[6,6] <- iratio('Os189Os192')[2]^2 # var(R92Os)
        E[7,7] <- iratio('Os190Os192')[2]^2 # var(R02Os)
        E[8,8] <- imass('Re')[2]^2 # var(MMRe)
        E[9,9] <- imass('Os')[2]^2 # var(MMOs)
    }
    nn <- dim(x$x)[1]
    out <- matrix(NA,nn,5) # X,sX,Y,sY,rXY with X=Re187/Os188 and Y=Os187/Os188
    colnames(out) <- c('X','sX','Y','sY','rXY')
    out[,1] <- Re187Os188
    out[,3] <- R78Os
    for (i in 1:nn){
        E[10,10] <- x$x[i,'errOs187Os188']^2 # var(R78Os) 
        E[11,11] <- (x$x[i,'errReppm'])^2    # var(Re)
        E[12,12] <- x$x[i,'errOsppm']^2      # var(Os)
        J[1,1] <- -Re187Os188[i]*dinvf187Re.dR57Re/invf187Re
        J[1,2] <- Re187Os188[i]*dinvf188Os.dR42Os/invf188Os[i]
        J[1,3] <- Re187Os188[i]*dinvf188Os.dR62Os/invf188Os[i]
        J[1,4] <- Re187Os188[i]*dinvf188Os.dR72Os/invf188Os[i]
        J[1,5] <- Re187Os188[i]*dinvf188Os.dR82Os/invf188Os[i]
        J[1,6] <- Re187Os188[i]*dinvf188Os.dR92Os/invf188Os[i]
        J[1,7] <- Re187Os188[i]*dinvf188Os.dR02Os/invf188Os[i]
        J[1,8] <- -Re187Os188[i]/MMRe
        J[1,9] <-  Re187Os188[i]/MMOs
        J[1,10] <- Re187Os188[i]*dinvf188Os.dR78Os/invf188Os[i]
        J[1,11] <-  Re187Os188[i]/Re[i]
        J[1,12] <- -Re187Os188[i]/Os[i]
        J[2,4] <-  dR78Os.dR72Os
        J[2,5] <-  dR78Os.dR82Os
        J[2,10] <- dR78Os.dR78Os
        E2 <- J %*% E %*% t(J)
        if (all(is.finite(E2))){
            out[i,2] <- sqrt(E2[1,1])
            out[i,4] <- sqrt(E2[2,2])
            out[i,5] <- stats::cov2cor(E2)[1,2]
        }
    }
    out
}

get.ReOs.ratio <- function(tt,st,exterr=TRUE){
    get.PD.ratio(tt,st,'Re187',exterr)
}

get.ReOs.age <- function(Os187Re187,sOs187Re187,exterr=TRUE){
    get.PD.age(Os187Re187,sOs187Re187,'Re187',exterr=exterr)
}

ReOs.age <- function(x,exterr=TRUE,i=NA,sigdig=NA,i2i=TRUE,...){
    PD.age(x,'Re187',exterr=exterr,i=i,sigdig=sigdig,i2i=i2i,...)
}
