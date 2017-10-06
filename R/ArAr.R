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
    Jac[1,1] <- (1-exp(L*tt))/J^2
    Jac[1,2] <- tt*exp(L*tt)/J
    Jac[1,3] <- L*exp(L*tt)/J
    E[1,1] <- sJ^2
    E[2,2] <- sL^2
    E[3,3] <- st^2
    sR <- sqrt(Jac %*% E %*% t(Jac))
    out <- c(R,sR)
}

# get a single Ar-Ar age
get.ArAr.age <- function(Ar40Ar39,sAr40Ar39,J,sJ,exterr=TRUE){
    L <- lambda("K40")[1]
    tt <- log(J*Ar40Ar39+1)/L
    Jac <- matrix(0,1,3)
    Jac[1,1] <- J/(L*(J*Ar40Ar39+1))
    Jac[1,2] <- Ar40Ar39/(L*(J*Ar40Ar39+1))
    if (exterr) Jac[1,3] <- -tt/L
    E <- matrix(0,3,3)
    E[1,1] <- sAr40Ar39^2
    E[2,2] <- sJ^2
    E[3,3] <- lambda("K40")[2]^2
    st <- sqrt(Jac %*% E %*% t(Jac))
    c(tt,st)
}

# x an object of class \code{ArAr} returns a matrix of 40Ar/39Ar-ages
# and their uncertainties. jcu = J-constant uncertainties.
ArAr.age <- function(x,jcu=TRUE,exterr=TRUE,i=NA,sigdig=NA,i2i=FALSE){
    ns <- length(x)
    if (ns<2) i2i <- FALSE
    out <- matrix(0,ns,2)
    colnames(out) <- c('t','s[t]')
    if (!jcu) x$J[2] <- 0
    rat <- ArAr.age.ratios(x)
    if (i2i){
        fit <- isochron.ArAr(x,plot=FALSE,exterr=exterr,inverse=FALSE)
        Ar4039x <- rat[,'Ar40Ar39'] - fit$a[1]/rat[,"Ar39Ar36"]
    } else {
        Ar4039x <- rat[,'Ar40Ar39'] - iratio("Ar40Ar36")[1]/rat[,"Ar39Ar36"]
    }
    J <- matrix(0,1,3)
    E <- matrix(0,3,3)
    for (j in 1:ns) {
        J[1,1] <- 1
        J[1,3] <- -1/rat[j,'Ar39Ar36']
        E[1,1] <- rat[j,'varAr40Ar39']
        E[2,2] <- rat[j,'varAr39Ar36']
        E[1,2] <- rat[j,'cov']
        E[2,1] <- E[1,2]
        if (i2i){
            J[1,2] <- fit$a[1]/rat[j,"Ar39Ar36"]^2
            if (exterr) E[3,3] <- fit$a[2]^2
        } else {
            J[1,2] <- iratio("Ar40Ar36")[1]/rat[j,"Ar39Ar36"]^2
            if (exterr) E[3,3] <- iratio("Ar40Ar36")[2]^2
        }
        sAr4039x <- sqrt(J %*% E %*% t(J))
        tt <- get.ArAr.age(Ar4039x[j],sAr4039x,x$J[1],x$J[2],exterr=exterr)
        out[j,] <- roundit(tt[1],tt[2],sigdig=sigdig)
    }
    if (!is.na(i)) out <- out[i,]
    out
}
