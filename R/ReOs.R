#' convert isotope dilution derived concentrations to ratios
#' x = matrix with columns 'Reppm','errReppm', 'Osppt','errOsppt' 
#' and 'Os187Os188','errOs187Os188'
#' @param exterr propagate the decay constant errors?
#' @param common remove the non-radiogenic component using the
#'     isotopic ratios stored in \code{settings('iratio',...)}
#' @noRd
ppm2ratios.ReOs <- function(x,exterr=FALSE,common=FALSE,...){
    R57Re <- iratio('Re185Re187')[1]
    R48Os <- iratio('Os184Os188')[1]
    R68Os <- iratio('Os186Os188')[1]
    R78Os0 <- iratio('Os187Os188')[1]
    R98Os <- iratio('Os189Os188')[1]
    R08Os <- iratio('Os190Os188')[1]
    R28Os <- iratio('Os192Os188')[1]
    if (common) { # make common Os correction
        R78Os <-  x$x[,'Os187Os188'] - R78Os0
        dR78Os.dR78Os0 <- -1
    } else {
        R78Os <- x$x[,'Os187Os188']
        dR78Os.dR78Os0 <- 0
    }
    Re <- x$x[,'Reppm']
    Os <- x$x[,'Osppm']
    MMOs <- imass('Os')[1]
    MMRe <- imass('Re')[1]
    invf187Re <- 1 + R57Re
    invf188Os <- R48Os + R68Os + R78Os + 1 + R98Os + R08Os + R28Os
    Re187Os188 <- (invf188Os/invf187Re)*(Re/Os)*(MMOs/MMRe)
    dinvf187Re.dR57Re <- 1
    dinvf188Os.dR48Os <- 1
    dinvf188Os.dR68Os <- 1
    dinvf188Os.dR78Os <- 1
    dinvf188Os.dR78Os0 <- dinvf188Os.dR78Os*dR78Os.dR78Os0
    dinvf188Os.dR98Os <- 1
    dinvf188Os.dR08Os <- 1
    dinvf188Os.dR28Os <- 1
    E <- matrix(0,12,12)
    J <- matrix(0,2,12)
    if (exterr){
        E[1,1] <- iratio('Re185Re187')[2]^2 # var(R57Re)
        E[2,2] <- iratio('Os184Os188')[2]^2 # var(R48Os)
        E[3,3] <- iratio('Os186Os188')[2]^2 # var(R68Os)
        E[4,4] <- iratio('Os187Os188')[2]^2 # var(R78Os0)
        E[5,5] <- iratio('Os189Os188')[2]^2 # var(R98Os)
        E[6,6] <- iratio('Os190Os188')[2]^2 # var(R08Os)
        E[7,7] <- iratio('Os192Os188')[2]^2 # var(R28Os)
        E[8,8] <- imass('Re')[2]^2
        E[9,9] <- imass('Os')[2]^2
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
        J[1,2] <-  Re187Os188[i]*dinvf188Os.dR48Os/invf188Os[i]
        J[1,3] <-  Re187Os188[i]*dinvf188Os.dR68Os/invf188Os[i]
        J[1,4] <-  Re187Os188[i]*dinvf188Os.dR78Os0/invf188Os[i]
        J[1,5] <-  Re187Os188[i]*dinvf188Os.dR98Os/invf188Os[i]
        J[1,6] <-  Re187Os188[i]*dinvf188Os.dR08Os/invf188Os[i]
        J[1,7] <-  Re187Os188[i]*dinvf188Os.dR28Os/invf188Os[i]
        J[1,8] <- -Re187Os188[i]/MMRe
        J[1,9] <-  Re187Os188[i]/MMOs
        J[1,10] <-  Re187Os188[i]*dinvf188Os.dR78Os/invf188Os[i]
        J[1,11] <-  Re187Os188[i]/Re[i]
        J[1,12] <- -Re187Os188[i]/Os[i]
        J[2,4] <-  dR78Os.dR78Os0
        J[2,10] <- 1
        E2 <- J %*% E %*% t(J)
        if (all(is.finite(E2))){
            out[i,2] <- sqrt(E2[1,1])
            out[i,4] <- sqrt(E2[2,2])
            out[i,5] <- stats::cov2cor(E2)[1,2]
        }
    }
    out
}

get_ReOs_ratio <- function(tt,st,exterr=FALSE){
    getDPratio(tt,st,'Re187',exterr)
}

get_ReOs_age <- function(Os187Re187,sOs187Re187,exterr=FALSE){
    getPDage(Os187Re187,sOs187Re187,'Re187',exterr=exterr)
}

ReOs_age <- function(x,exterr=FALSE,i=NULL,i2i=TRUE,projerr=FALSE,...){
    PD_age(x,'Re187',exterr=exterr,i=i,i2i=i2i,projerr=projerr,...)
}
