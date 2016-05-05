get.covmat.UPb <- function(X,i){
    covmat <- matrix(rep(0,9),nrow=3)
    if (X$format == 1){
        rownames(covmat) <- c('Pb207Pb206','Pb206U238','Pb207U235')
        colnames(covmat) <- c('Pb207Pb206','Pb206U238','Pb207U235')
        relvar207 <- 0.5 * ((X$x[i,'errPb207Pb206']/X$x[i,'Pb207Pb206'])^2 +
                            (X$x[i,'errPb207U235']/X$x[i,'Pb207U235'])^2 -
                            (X$x[i,'errPb206U238']/X$x[i,'Pb206U238'])^2)
        relvar206 <- 0.5 * ((X$x[i,'errPb207Pb206']/X$x[i,'Pb207Pb206'])^2 -
                            (X$x[i,'errPb207U235']/X$x[i,'Pb207U235'])^2 +
                            (X$x[i,'errPb206U238']/X$x[i,'Pb206U238'])^2)
        relvar238 <- 0.5 * ((X$x[i,'errPb207U235']/X$x[i,'Pb207U235'])^2 +
                            (X$x[i,'errPb206U238']/X$x[i,'Pb206U238'])^2 -
                            (X$x[i,'errPb207Pb206']/X$x[i,'Pb207Pb206'])^2)
        covmat[1,1] <- X$x[i,'errPb207Pb206']^2
        covmat[2,2] <- X$x[i,'errPb206U238']^2
        covmat[3,3] <- X$x[i,'errPb207U235']^2
        covmat[1,2] <- -relvar206*X$x[i,'Pb207Pb206']*X$x[i,'Pb206U238']
        covmat[1,3] <- relvar207/(X$x[i,'Pb207Pb206']*X$x[i,'Pb207U235'])
        covmat[2,3] <- relvar238*X$x[i,'Pb206U238']*X$x[i,'Pb207U235']
        covmat[2,1] <- covmat[1,2]
        covmat[3,2] <- covmat[2,3]
        covmat[3,1] <- covmat[1,3]
    }
    colnames(covmat) <- rownames(covmat)
    covmat
}

get.ratios.UPb <- function(age){
    if (age == 0){ age <- 1e-10 }
    out <- list()
    l8 <- lambda('U238')$x
    l5 <- lambda('U235')$x
    R.x <- 1/R238235()$x
    R.e <- R.x*R238235()$e/R238235()$x
    Pb206U238 <- (exp(l8*age)-1)
    U238Pb206 <- 1/(exp(l8*age)-1)
    Pb207U235 <- (exp(l5*age)-1)
    Pb207Pb206 <- R.x*Pb207U235/Pb206U238
    out$x <- c(Pb207Pb206,Pb206U238,Pb207U235,U238Pb206)
    E <- matrix(0,nrow=3,ncol=3)
    E[1,1] <- R.e^2
    E[2,2] <- lambda('U238')$e^2
    E[3,3] <- lambda('U235')$e^2
    J <- matrix(0,nrow=4,ncol=3)
    J[1,1] <- (exp(l5*age)-1)/(exp(l8*age)-1)
    J[1,2] <- R.x*age*exp(l5*age)/(exp(l8*age)-1)
    J[1,3] <- -R.x*(exp(l5*age)-1)*age*exp(l8*age)/(exp(l8*age)-1)^2
    J[2,2] <- age*exp(l8*age)
    J[3,3] <- age*exp(l5*age)
    J[4,2] <- age*exp(l8*age)/(exp(l8*age)-1)^2
    out$cov <- J %*% E %*% t(J)
    names(out$x) <- c('Pb207Pb206','Pb206U238','Pb207U235','U238Pb206')
    rownames(out$cov) <- names(out$x)
    colnames(out$cov) <- names(out$x)
    out
}

get.ages.UPb <- function(x){
    
}

# returns a five item list containing \code{x}: the geometric mean
# isotopic composition, \code{cov}: the covariance matrix of the
# geometric mean isotopic composition, \code{mswd}: the reduced
# chi-square statistic for the geometric mean composition,
# \code{age}: the 206Pb/238U-207Pb/235U concordia age, and
# \code{err}: its standard error
concordia.age <- function(x){

}

discordia.age <- function(x){

}

get.Pb207U235age <- function(Pb207U235){
    log(1+Pb207U235)/lambda('U235')$x
}

get.Pb206U238age <- function(Pb206U238){
    log(1+Pb206U238)/lambda('U238')$x
}

get.Pb207Pb206age <- function(Pb207Pb206){
    Pb207Pb206.misfit <- function(x,y) { (get.ratios.UPb(x)$x[1] - y)^2 }
    out <- stats::optimize(Pb207Pb206.misfit,c(0,4600),y=Pb207Pb206)
    out$minimum
}
