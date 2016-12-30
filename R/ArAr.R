get.selection.ArAr <- function(x,inverse=TRUE,...){
    if (inverse) selection <- c('Ar39Ar40','Ar36Ar40')
    else selection <- c('Ar39Ar36','Ar40Ar36')
    selection
}

get.covmat.ArAr <- function(x,i,...){
    out <- matrix(rep(0,16),nrow=4)
    rownames(out) <- c('Ar39Ar40','Ar36Ar40','Ar39Ar36','Ar40Ar36')
    colnames(out) <- rownames(out)
    if (x$format == 1 | x$format == 2){
        out['Ar39Ar40','Ar39Ar40'] <- x$x[i,'errAr39Ar40']^2
        out['Ar36Ar40','Ar36Ar40'] <- x$x[i,'errAr36Ar40']^2
        out['Ar39Ar36','Ar39Ar36'] <- x$x[i,'errAr39Ar36']^2
        out['Ar40Ar36','Ar40Ar36'] <- x$x[i,'errAr40Ar36']^2
        
        out['Ar39Ar40','Ar36Ar40'] <- get.cov.xzyz(
            x$x[i,'Ar39Ar40'],x$x[i,'errAr39Ar40'],
            x$x[i,'Ar36Ar40'],x$x[i,'errAr36Ar40'],x$x[i,'errAr39Ar36']
        )
        out['Ar39Ar40','Ar39Ar36'] <- get.cov.zxzy(
            x$x[i,'Ar39Ar40'],x$x[i,'errAr39Ar40'],
            x$x[i,'Ar39Ar36'],x$x[i,'errAr39Ar36'],x$x[i,'errAr40Ar36']
        )
        out['Ar39Ar40','Ar40Ar36'] <- get.cov.xzzy(
            x$x[i,'Ar39Ar40'],x$x[i,'errAr39Ar40'],
            x$x[i,'Ar40Ar36'],x$x[i,'errAr40Ar36'],x$x[i,'errAr39Ar36']
        )
        out['Ar36Ar40','Ar39Ar36'] <- get.cov.xzzy(
            x$x[i,'Ar39Ar36'],x$x[i,'errAr39Ar36'],
            x$x[i,'Ar36Ar40'],x$x[i,'errAr36Ar40'],x$x[i,'errAr39Ar40']
        )
        out['Ar36Ar40','Ar40Ar36'] <- 0
        out['Ar39Ar36','Ar40Ar36'] <- get.cov.xzyz(
            x$x[i,'Ar39Ar36'],x$x[i,'errAr39Ar36'],
            x$x[i,'Ar40Ar36'],x$x[i,'errAr40Ar36'],x$x[i,'errAr39Ar40']
        )
        
        out['Ar36Ar40','Ar39Ar40'] <- out['Ar39Ar40','Ar36Ar40']
        out['Ar39Ar36','Ar39Ar40'] <- out['Ar39Ar40','Ar39Ar36']
        out['Ar40Ar36','Ar39Ar40'] <- out['Ar39Ar40','Ar40Ar36']
        out['Ar39Ar36','Ar36Ar40'] <- out['Ar36Ar40','Ar39Ar36']
        out['Ar40Ar36','Ar36Ar40'] <- out['Ar36Ar40','Ar40Ar36']
        out['Ar40Ar36','Ar39Ar36'] <- out['Ar39Ar36','Ar40Ar36']
    }
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
    if (!jcu) x$J[2] <- 0
    ns <- nrow(x$x)
    out <- matrix(0,ns,2)
    colnames(out) <- c('t','s[t]')
    Ar4036 <- x$x[,"Ar40Ar36"]
    sAr4036 <- x$x[,"errAr40Ar36"]
    Ar4039 <- 1/x$x[,"Ar39Ar40"]
    sAr4039 <- x$x[,"errAr39Ar40"]/x$x[,"Ar39Ar40"]^2
    if (i2i){
        fit <- isochron.ArAr(x,plot=FALSE,exterr=exterr,inverse=FALSE)
        Ar4039x <- Ar4039 - fit$a[1]/x$x[,"Ar39Ar36"]
    } else {
        Ar4039x <- Ar4039-iratio("Ar40Ar36")[1]/x$x[,"Ar39Ar36"]
    }
    J <- matrix(0,1,3)
    E <- matrix(0,3,3)
    for (j in 1:ns) {
        covmat <- get.covmat.ArAr(x,j)
        J[1,1] <- -1/x$x[j,"Ar39Ar40"]^2
        J[1,3] <- -1/x$x[j,"Ar39Ar36"]
        E[1:2,1:2] <- covmat[c("Ar39Ar40","Ar39Ar36"),c("Ar39Ar40","Ar39Ar36")]
        if (i2i) {
            J[1,2] <- fit$a[1]/x$x[j,"Ar39Ar36"]^2
            if (exterr) E[3,3] <- fit$a[2]^2
        } else {
            J[1,2] <- iratio("Ar40Ar36")[1]/x$x[j,"Ar39Ar36"]^2
            if (exterr) E[3,3] <- iratio("Ar40Ar36")[2]^2
        }
        sAr4039x <- sqrt(J %*% E %*% t(J))
        tt <- get.ArAr.age(Ar4039x[j],sAr4039x,x$J[1],x$J[2],exterr=exterr)
        t.out <- roundit(tt[1],tt[2],sigdig=sigdig)
        out[j,] <- c(t.out$x,t.out$err)
    }
    if (!is.na(i)) out <- out[i,]
    out
}
