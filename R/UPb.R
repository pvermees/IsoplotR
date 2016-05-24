get.covmat.UPb <- function(x,i){
    covmat <- matrix(rep(0,16),nrow=4)
    rownames(covmat) <- c('Pb207U235','Pb206U238','U238Pb206','Pb207Pb206')
    colnames(covmat) <- rownames(covmat)
    if (x$format == 1){

        covmat[1,1] <- x$x[i,'errPb207U235']^2
        covmat[2,2] <- x$x[i,'errPb206U238']^2
        covmat[3,3] <- (x$x[i,'errPb206U238']^2)/(x$x[i,'Pb206U238']^4)
        covmat[4,4] <- x$x[i,'errPb207Pb206']^2

        covmat[1,2] <- 0.5 * x$x[i,'Pb207U235'] * x$x[i,'Pb206U238'] * (
                            (x$x[i,'errPb207U235']/x$x[i,'Pb207U235'])^2 +
                            (x$x[i,'errPb206U238']/x$x[i,'Pb206U238'])^2 -
                            (x$x[i,'errPb207Pb206']/x$x[i,'Pb207Pb206'])^2 )
        covmat[1,3] <- 0.5 * x$x[i,'Pb207U235'] * x$x[i,'Pb206U238'] * (
                            (x$x[i,'errPb207Pb206']/x$x[i,'Pb207Pb206'])^2 -
                            (x$x[i,'errPb207U235']/x$x[i,'Pb207U235'])^2 -
                            covmat[3,3]*x$x[i,'Pb206U238']^2 )
        covmat[1,4] <- 0.5 * x$x[i,'Pb207U235']*x$x[i,'Pb207Pb206'] * (
                            (x$x[i,'errPb207U235']/x$x[i,'Pb207U235'])^2 +
                            (x$x[i,'errPb207Pb206']/x$x[i,'Pb207Pb206'])^2 -
                            (x$x[i,'errPb206U238']/x$x[i,'Pb206U238'])^2 )
        covmat[2,1] <- covmat[1,2]
        covmat[2,3] <- -1
        covmat[2,4] <- 0.5 * x$x[i,'Pb207Pb206']*x$x[i,'Pb206U238'] * (
                            (x$x[i,'errPb207U235']/x$x[i,'Pb207U235'])^2 -
                            (x$x[i,'errPb207Pb206']/x$x[i,'Pb207Pb206'])^2 -
                            (x$x[i,'errPb206U238']/x$x[i,'Pb206U238'])^2 )
        covmat[3,1] <- covmat[1,3]
        covmat[3,2] <- -1
        covmat[3,4] <- 0.5 * x$x[i,'Pb207Pb206']/x$x[i,'Pb206U238'] * (
                            (x$x[i,'errPb207Pb206']/x$x[i,'Pb207Pb206'])^2 +
                             covmat[3,3]*x$x[i,'Pb206U238']^2 -
                            (x$x[i,'errPb207U235']/x$x[i,'Pb207U235'])^2 )
        covmat[4,1] <- covmat[1,4]
        covmat[4,2] <- covmat[2,4]
        covmat[4,3] <- covmat[3,4]
    }
    covmat
}

get.ratios.UPb <- function(age){
    if (age == 0){ age <- 1e-10 }
    out <- list()
    l8 <- lambda('U238')[1]
    l5 <- lambda('U235')[1]
    R.x <- 1/I.R('U238U235')[1]
    R.e <- R.x*I.R('U238U235')[2]/I.R('U238U235')[1]

    Pb207U235 <- (exp(l5*age)-1)
    Pb206U238 <- (exp(l8*age)-1)
    U238Pb206 <- 1/(exp(l8*age)-1)
    Pb207Pb206 <- R.x*Pb207U235/Pb206U238
    
    out$x <- c(Pb207U235,Pb206U238,U238Pb206,Pb207Pb206)

    E <- matrix(0,3,3)
    E[1,1] <- R.e^2
    E[2,2] <- lambda('U235')[2]^2
    E[3,3] <- lambda('U238')[2]^2

    J <- matrix(0,nrow=4,ncol=3)
    J[1,2] <- age*exp(l5*age)
    J[2,3] <- age*exp(l8*age)
    J[3,3] <- -age*exp(l8*age)/(exp(l8*age)-1)^2
    J[4,1] <- (exp(l5*age)-1)/(exp(l8*age)-1)
    J[4,2] <- R.x*age*exp(l5*age)/(exp(l8*age)-1)
    J[4,3] <- -R.x*(exp(l5*age)-1)*age*exp(l8*age)/(exp(l8*age)-1)^2
    out$cov <- J %*% E %*% t(J)
    names(out$x) <- c('Pb207U235','Pb206U238','U238Pb206','Pb207Pb206')
    rownames(out$cov) <- names(out$x)
    colnames(out$cov) <- names(out$x)
    out
}

get.Pb207U235age <- function(Pb207U235){
    log(1+Pb207U235)/lambda('U235')[1]
}

get.Pb206U238age <- function(Pb206U238){
    log(1+Pb206U238)/lambda('U238')[1]
}

#' @importFrom stats optimize
get.Pb207Pb206age <- function(Pb207Pb206){
    Pb207Pb206.misfit <- function(x,y) { (get.ratios.UPb(x)$x[1] - y)^2 }
    out <- optimize(Pb207Pb206.misfit,c(0,4600),y=Pb207Pb206)
    out$minimum
}

UPb.preprocess <- function(x,wetherill){
    selection <- get.UPb.selection(wetherill)
    out <- list()
    for (i in 1:nrow(x$x)){
        X <- x$x[i,selection]
        covmat <- get.covmat.UPb(x,i)[selection,selection]
        out[[i]] <- list(x=X, cov=covmat)
    }
    out   
}
