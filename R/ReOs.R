# convert isotope dilution derived concentrations to ratios
# x = matrix with columns 'Reppm','errReppm', 'Osppt','errOsppt' 
# and 'Os187Os188','errOs187Os188'
ID.Re <- function(x,exterr=FALSE,isochron=TRUE){
    R57Re <- iratio('Re185Re187')[1]
    R42Os <- iratio('Os184Os192')[1]
    R62Os <- iratio('Os186Os192')[1]
    R72Os <- iratio('Os187Os192')[1]
    R82Os <- iratio('Os188Os192')[1]
    R02Os <- iratio('Os190Os192')[1]
    R92Os <- iratio('Os189Os192')[1]
    dR78Os.dR78Os <- 1
    if (isochron) {
        R78Os <- x$x[,'Os187Os188']
        dR78Os.dR72Os <- 0
        dR78Os.dR82Os <- 0
    } else { # make common Os correction
        R78Os <-  x$x[,'Os187Os188'] - R72Os/R82Os
        dR78Os.dR72Os <- -1/R82Os
        dR78Os.dR82Os <- R72Os/R82Os^2 
    }
    Re <- x$x[,'Reppm']*1000
    Os <- x$x[,'Osppt']
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
        E[8,8] <- imass('Re')[2]^2
        E[9,9] <- imass('Os')[2]^2
    }
    nn <- dim(x$x)[1]
    out <- matrix(NA,nn,5) # X,sX,Y,sY,rXY with X=Re187/Os188 and Y=Os187/Os188
    colnames(out) <- c('X','sX','Y','sY','rXY')
    out[,1] <- Re187Os188
    out[,3] <- R78Os
    for (i in 1:nn){
        E[10,10] <- x$x[i,'errOs187Os188']^2
        E[11,11] <- (x$x[i,'errReppm']*1000)^2
        E[12,12] <- x$x[i,'errOsppt']^2
        J[1,1] <- (df187Re.dR57Re*invf188Os[i])*(Re[i]/Os[i])*(MMOs/MMRe)
        J[1,2] <- (f187Re*dinvf188Os.dR42Os)*(Re[i]/Os[i])*(MMOs/MMRe)
        J[1,3] <- (f187Re*dinvf188Os.dR62Os)*(Re[i]/Os[i])*(MMOs/MMRe)
        J[1,5] <- (f187Re*dinvf188Os.dR82Os)*(Re[i]/Os[i])*(MMOs/MMRe)
        J[1,6] <- (f187Re*dinvf188Os.dR92Os)*(Re[i]/Os[i])*(MMOs/MMRe)
        J[1,7] <- (f187Re*dinvf188Os.dR02Os)*(Re[i]/Os[i])*(MMOs/MMRe)
        J[1,8] <- -Re187Os188[i]/MMRe
        J[1,9] <- (f187Re*invf188Os[i])*(Re[i]/Os[i])*(1/MMRe)
        J[1,10] <- (f187Re*dinvf188Os.dR78Os)*(Re[i]/Os[i])*(MMOs/MMRe)
        J[1,11] <- (f187Re*invf188Os[i])*(1/Os[i])*(MMOs/MMRe)
        J[1,12] <- -Re187Os188[i]/Re[i]
        J[2,4] <- dR78Os.dR72Os
        J[2,5] <- dR78Os.dR82Os
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

get.ReOs.age <- function(Os187Re187,sOs187Re187,exterr=TRUE){
    l187 <- lambda('Re187')
    tt <- log(1 + Os187Re187)/l187[1]
    E <- matrix(0,2,2)
    J <- matrix(0,1,2)
    E[1,1] <- sOs187Re187^2
    if (exterr) E[2,2] <- l187[2]^2
    J[1,1] <- 1/(l187[1]*(1 + Os187Re187))
    J[1,2] <- -tt/l187[1]
    st <- sqrt(J %*% E %*% t(J))
    c(tt,st)
}

ReOs.age <- function(x,jcu=TRUE,exterr=TRUE,i=NA,sigdig=NA){
    ns <- nrow(x$x)
    dat <- ID.Re(x,exterr=exterr,isochron=FALSE)
    out <- matrix(0,ns,2)
    E <- matrix(0,2,2)
    J <- matrix(0,1,2)
    for (j in 1:ns) {
        E[1,1] <- dat[j,'sX']^2
        E[2,2] <- dat[j,'sY']^2
        E[1,2] <- dat[j,'rXY']*dat[j,'sX']*dat[j,'sY']
        E[2,1] <- E[1,2]
        J[1,1] <- 1/dat[j,'Y']
        J[1,2] <- -dat[j,'X']/dat[j,'Y']^2        
        Os187Re187 <- dat[j,'Y']/dat[j,'X']
        sOs187Re187 <- sqrt(J%*%E%*%t(J))
        out[j,] <- get.ReOs.age(Os187Re187,sOs187Re187,exterr=TRUE)
    }
    if (!is.na(i)) out <- out[i,]
    out
}
