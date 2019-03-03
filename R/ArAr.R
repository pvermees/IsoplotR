# returns ratios, errors and correlations
ArAr.normal.ratios <- function(x){
    ns <- length(x)
    out <- matrix(0,ns,5)
    if (x$format==1){
        Ar39Ar36 <- x$x[,'Ar39Ar36']
        errAr39Ar36 <- x$x[,'errAr39Ar36']
        Ar40Ar36 <- x$x[,'Ar40Ar36']
        errAr40Ar36 <- x$x[,'errAr40Ar36']
        rho <- x$x[,'rho']
    } else if (x$format==2){
        Ar39Ar36 <- x$x[,'Ar39Ar40']/x$x[,'Ar36Ar40']
        Ar40Ar36 <- 1/x$x[,'Ar36Ar40']
        J11 <- 1/x$x[,'Ar36Ar40']
        J22 <- -1/x$x[,'Ar36Ar40']^2
        J12 <- -x$x[,'Ar39Ar40']/x$x[,'Ar36Ar40']^2
        J21 <- 0
        E11 <- x$x[,'errAr39Ar40']^2
        E22 <- x$x[,'errAr36Ar40']^2
        E12 <- x$x[,'rho']*x$x[,'errAr39Ar40']*x$x[,'errAr36Ar40']
        covmat <- errorprop(J11,J12,J21,J22,E11,E22,E12)
        errAr39Ar36 <- sqrt(covmat[,'varX'])
        errAr40Ar36 <- sqrt(covmat[,'varY'])
        rho <- covmat[,'cov']/(errAr39Ar36*errAr40Ar36)
    } else if (x$format==3){
        Ar39Ar36 <- x$x[,'Ar39Ar36']
        Ar40Ar36 <- 1/x$x[,'Ar36Ar40']
        errAr39Ar36 <- x$x[,'errAr39Ar36']
        errAr40Ar36 <- x$x[,'errAr36Ar40']/x$x[,'Ar36Ar40']^2
        rho <- get.cor.div(Ar39Ar36,errAr39Ar36,
                           Ar40Ar36,errAr40Ar36,
                           x$x[,'Ar39Ar40'],x$x[,'errAr39Ar40'])
    }
    out <- cbind(Ar39Ar36,errAr39Ar36,Ar40Ar36,errAr40Ar36,rho)
    colnames(out) <- c('Ar39Ar36','errAr39Ar36',
                       'Ar40Ar36','errAr40Ar36','rho')
    out
}
# returns ratios, errors and correlations
ArAr.inverse.ratios <- function(x){
    ns <- length(x)
    out <- matrix(0,ns,5)
    if (x$format==1){
        Ar39Ar40 <- x$x[,'Ar39Ar36']/x$x[,'Ar40Ar36']
        Ar36Ar40 <- 1/x$x[,'Ar40Ar36']
        J11 <- 1/x$x[,'Ar40Ar36']
        J22 <- -1/x$x[,'Ar40Ar36']^2
        J12 <- -x$x[,'Ar39Ar36']/x$x[,'Ar40Ar36']^2
        J21 <- 0
        E11 <- x$x[,'errAr39Ar36']^2
        E22 <- x$x[,'errAr40Ar36']^2
        E12 <- x$x[,'rho']*x$x[,'errAr39Ar36']*x$x[,'errAr40Ar36']
        covmat <- errorprop(J11,J12,J21,J22,E11,E22,E12)
        errAr39Ar40 <- sqrt(covmat[,'varX'])
        errAr36Ar40 <- sqrt(covmat[,'varY'])
        rho <- covmat[,'cov']/(errAr39Ar40*errAr36Ar40)
    } else if (x$format==2){
        Ar39Ar40 <- x$x[,'Ar39Ar40']
        errAr39Ar40 <- x$x[,'errAr39Ar40']
        Ar36Ar40 <- x$x[,'Ar36Ar40']
        errAr36Ar40 <- x$x[,'errAr36Ar40']
        rho <- x$x[,'rho']
    } else if (x$format==3){
        Ar39Ar40 <- x$x[,'Ar39Ar40']
        Ar36Ar40 <- x$x[,'Ar36Ar40']
        errAr39Ar40 <- x$x[,'errAr39Ar40']
        errAr36Ar40 <- x$x[,'errAr36Ar40']
        rho <- get.cor.div(Ar39Ar40,errAr39Ar40,
                           Ar36Ar40,errAr36Ar40,
                           x$x[,'Ar39Ar36'],x$x[,'errAr39Ar36'])
    }
    out <- cbind(Ar39Ar40,errAr39Ar40,Ar36Ar40,errAr36Ar40,rho)
    colnames(out) <- c('Ar39Ar40','errAr39Ar40',
                       'Ar36Ar40','errAr36Ar40','rho')
    out
}
# returns ratios, variances and covariances
ArAr.age.ratios <- function(x){
    ns <- length(x)
    out <- matrix(0,ns,5)
    if (x$format==1){
        Ar40Ar39 <- x$x[,'Ar40Ar36']/x$x[,'Ar39Ar36']
        Ar39Ar36 <- x$x[,'Ar39Ar36']
        J11 <- -x$x[,'Ar40Ar36']/x$x[,'Ar39Ar36']^2
        J22 <- 0
        J12 <- 1/x$x[,'Ar39Ar36']
        J21 <- 1
        E11 <- x$x[,'errAr39Ar36']^2
        E22 <- x$x[,'errAr40Ar36']^2
        E12 <- x$x[,'rho']*x$x[,'errAr39Ar36']*x$x[,'errAr40Ar36']
        covmat <- errorprop(J11,J12,J21,J22,E11,E22,E12)
    } else if (x$format==2){
        Ar40Ar39 <- 1/x$x[,'Ar39Ar40']
        Ar39Ar36 <- x$x[,'Ar39Ar40']/x$x[,'Ar36Ar40']
        J11 <- -1/x$x[,'Ar39Ar40']^2
        J22 <- -x$x[,'Ar39Ar40']/x$x[,'Ar36Ar40']^2
        J12 <- 0
        J21 <- 1/x$x[,'Ar36Ar40']
        E11 <- x$x[,'errAr39Ar40']^2
        E22 <- x$x[,'errAr36Ar40']^2
        E12 <- x$x[,'rho']*x$x[,'errAr39Ar40']*x$x[,'errAr36Ar40']
        covmat <- errorprop(J11,J12,J21,J22,E11,E22,E12)
    } else if (x$format==3){
        Ar40Ar39 <- 1/x$x[,'Ar39Ar40']
        Ar39Ar36 <- x$x[,'Ar39Ar36']
        errAr40Ar39 <- x$x[,'errAr39Ar40']/x$x[,'Ar39Ar40']^2
        errAr39Ar36 <- x$x[,'errAr39Ar36']
        errAr40Ar36 <- x$x[,'errAr36Ar40']/x$x[,'Ar36Ar40']^2
        covmat <- matrix(0,ns,3)
        colnames(covmat) <- c('varX','varY','cov')
        covmat[,'varX'] <- errAr40Ar39^2
        covmat[,'varY'] <- errAr39Ar36^2
        covmat[,'cov'] <- get.cov.mult(Ar40Ar39,errAr40Ar39,
                                       Ar39Ar36,errAr39Ar36,
                                       Ar40Ar39*Ar39Ar36,errAr40Ar36)
    }
    out <- cbind(Ar40Ar39,covmat[,'varX'],Ar39Ar36,covmat[,'varY'],covmat[,'cov'])
    colnames(out) <- c('Ar40Ar39','varAr40Ar39',
                       'Ar39Ar36','varAr39Ar36','cov')
    out
}

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
ArAr.age <- function(x,exterr=TRUE,i=NA,sigdig=NA,i2i=FALSE){
    ns <- length(x)
    if (ns<2) i2i <- FALSE
    out <- matrix(0,ns,2)
    colnames(out) <- c('t','s[t]')
    y <- data2york(x,inverse=TRUE)
    if (i2i){
        fit <- york(y)
        b <- fit$b[1]
    } else {
        b <- iratio("Ar40Ar36")[1]
    }
    DP <- 1/(y[,'X'] - y[,'Y']/b)
    J1 <- -(DP^2)
    J2 <- (DP^2)/b
    E11 <- y[,'sX']^2
    E22 <- y[,'sY']^2
    E12 <- y[,'rXY']*y[,'sX']*y[,'sY']
    if (exterr){
        J3 <- -(DP^2)*y[,'sY']/b^2
        E33 <- iratio("Ar40Ar36")[2]^2
        sDP <- errorprop1x3(J1,J2,J3,E11,E22,E33,E12)
    } else {
        sDP <- errorprop1x2(J1,J2,E11,E22,E12)
    }
    tt <- get.ArAr.age(DP,sDP,x$J[1],x$J[2],exterr=exterr)
    out <- roundit(tt[,1],tt[,2],sigdig=sigdig)
    if (!is.na(i)) out <- out[i,]
    out
}
