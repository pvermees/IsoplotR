# x = matrix with columns 'Reppm','errReppm', 'Osppt','errOsppt' 
# and 'Os187Os188','errOs187Os188')
ID.Re <- function(x){
    R57Re <- iratio('Re185Re187')[1]
    R42Os <- iratio('Os184Os192')[1]
    R62Os <- iratio('Os186Os192')[1]
    R82Os <- iratio('Os188Os192')[1]
    R02Os <- iratio('Os190Os192')[1]
    R92Os <- iratio('Os189Os192')[1]
    R78Os <- x[,'Os187Os188']
    Re <- x[,'Reppm']*1000
    Os <- x[,'Osppt']
    MMOs <- imass('Os')[1]
    MMRe <- imass('Re')[1]
    invf188Os <- R42Os/R82Os + R62Os/R82Os + R78Os + 1 +
                 R92Os/R82Os + R02Os/R82Os + 1/R82Os
    f187Re <- 1/(1+R57Re)
    Re187Os188 <- (f187Re*invf188Os)*(Re/Os)*(MMOs/MMRe)
    df187Re.dR57Re <- -1/(1+R57Re)^2
    dinvf188Os.dR42Os <- 1/R82Os
    dinvf188Os.dR62Os <- 1/R82Os
    dinvf188Os.dR78Os <- 1
    dinvf188Os.dR82Os <- -(R42Os + R82Os + R92Os + 1)/R82Os^2
    dinvf188Os.dR02Os <- 1/R82Os
    dinvf188Os.dR92Os <- 1/R82Os
    E <- matrix(0,11,11)
    J <- matrix(0,2,11)
    E[1,1] <- iratio('Re185Re187')[2]^2 # var(R57Re)
    E[2,2] <- iratio('Os184Os192')[2]^2 # var(R42Os)
    E[3,3] <- iratio('Os186Os192')[2]^2 # var(R62Os)
    E[4,4] <- iratio('Os188Os192')[2]^2 # var(R82Os)
    E[5,5] <- iratio('Os189Os192')[2]^2 # var(R92Os)
    E[6,6] <- iratio('Os190Os192')[2]^2 # var(R02Os)
    E[7,7] <- imass('Re')[2]^2
    E[8,8] <- imass('Os')[2]^2
    nn <- dim(x)[1]
    out <- matrix(0,nn,5) # X,sX,Y,sY,rXY with X=Re187/Os188 and Y=Os187/Os188
    colnames(out) <- c('X','sX','Y','sY','rXY')
    out[,1] <- Re187Os188
    out[,3] <- R78Os
    out[,4] <- x[,'errOs187Os188']
    for (i in 1:nn){
        E[9,9] <- x[i,'errOs187Os188']^2
        E[10,10] <- (x[i,'errReppm']*1000)^2
        E[11,11] <- x[i,'errOsppt']^2
        J[1,1] <- (df187Re.dR57Re*invf188Os[i])*(Re[i]/Os[i])*(MMOs/MMRe)
        J[1,2] <- (f187Re*dinvf188Os.dR42Os)*(Re[i]/Os[i])*(MMOs/MMRe)
        J[1,3] <- (f187Re*dinvf188Os.dR62Os)*(Re[i]/Os[i])*(MMOs/MMRe)
        J[1,4] <- (f187Re*dinvf188Os.dR82Os)*(Re[i]/Os[i])*(MMOs/MMRe)
        J[1,5] <- (f187Re*dinvf188Os.dR92Os)*(Re[i]/Os[i])*(MMOs/MMRe)
        J[1,6] <- (f187Re*dinvf188Os.dR02Os)*(Re[i]/Os[i])*(MMOs/MMRe)
        J[1,7] <- -Re187Os188[i]/MMRe
        J[1,8] <- (f187Re*invf188Os[i])*(Re[i]/Os[i])*(1/MMRe)
        J[1,9] <- (f187Re*dinvf188Os.dR78Os)*(Re[i]/Os[i])*(MMOs/MMRe)
        J[1,10] <- (f187Re*invf188Os[i])*(1/Os[i])*(MMOs/MMRe)
        J[1,11] <- -Re187Os188[i]/Re[i]
        J[2,9] <- 1
        E2 <- J %*% E %*% t(J)
        out[i,2] <- sqrt(E2[1,1])
        out[i,5] <- stats::cov2cor(E2)[1,2]
    }
    out
}

get.ReOs.Age <- function(x){
    
}
