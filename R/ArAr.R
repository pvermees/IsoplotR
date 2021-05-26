get.ArAr.ratio <- function(tt,st,J,sJ,exterr=TRUE){
    L <- lambda("K40")[1]
    sL <- lambda("K40")[2]
    R <- (exp(L*tt)-1)/J
    Jac <- matrix(0,1,3)
    E <- matrix(0,3,3)
    if (exterr){
        Jac[1,1] <- (1-exp(L*tt))/J^2
        Jac[1,2] <- tt*exp(L*tt)/J
    }
    Jac[1,3] <- L*exp(L*tt)/J
    E[1,1] <- sJ^2
    E[2,2] <- sL^2
    E[3,3] <- st^2
    sR <- sqrt(Jac %*% E %*% t(Jac))
    out <- c(R,sR)
}

get.ArAr.age <- function(Ar40Ar39,sAr40Ar39,J,sJ,exterr=TRUE){
    L <- lambda("K40")[1]
    tt <- log(J*Ar40Ar39+1)/L
    J1 <- J/(L*(J*Ar40Ar39+1))
    if (exterr){
        J2 <- Ar40Ar39/(L*(J*Ar40Ar39+1))
        J3 <- -tt/L
    } else {
        J2 <- 0
        J3 <- 0
    }
    E11 <- sAr40Ar39^2
    E22 <- sJ^2
    E33 <- lambda("K40")[2]^2
    st <- errorprop1x3(J1,J2,J3,E11,E22,E33)
    out <- cbind(tt,st)
    colnames(out) <- c('t','st')
    out
}

# x an object of class \code{ArAr} returns a matrix 
# of 40Ar/39Ar-ages and their uncertainties.
ArAr.age <- function(x,exterr=TRUE,i=NA,sigdig=NA,
                     i2i=FALSE,projerr=FALSE,omit4c=NULL){
    ns <- length(x)
    if (ns<2) i2i <- FALSE
    out <- matrix(0,ns,2)
    colnames(out) <- c('t','s[t]')
    if (i2i){
        y <- data2york(x,inverse=TRUE)
        fit <- york(clear(y,omit4c))
        b <- fit$b
        DP <- 1/(y[,'X'] - y[,'Y']/b[1])
        J1 <- -DP^2
        J2 <- (DP^2)/b[1]
        J3 <- -y[,'Y']*(DP/b[1])^2
        E33 <- rep(b[2]^2,ns)
    } else {
        y <- data2york(x,inverse=FALSE)
        a <- iratio("Ar40Ar36")
        DP <- (y[,'Y']-a[1])/y[,'X']
        J1 <- -DP/y[,'X']
        J2 <- 1/y[,'X']
        J3 <- -1/y[,'X']
        E33 <- rep(a[2]^2,ns)
    }
    E11 <- y[,'sX']^2
    E22 <- y[,'sY']^2
    E12 <- y[,'rXY']*y[,'sX']*y[,'sY']
    if (projerr) sDP <- errorprop1x3(J1,J2,J3,E11,E22,E33,E12)
    else sDP <- errorprop1x2(J1,J2,E11,E22,E12)
    tt <- get.ArAr.age(DP,sDP,x$J[1],x$J[2],exterr=exterr)
    out <- roundit(tt[,1],tt[,2],sigdig=sigdig)
    if (!is.na(i)) out <- out[i,]
    out
}
