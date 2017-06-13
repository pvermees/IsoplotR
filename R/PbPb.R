# returns ratios, errors and correlations
PbPb.normal.ratios <- function(x){
    ns <- length(x)
    out <- matrix(0,ns,5)
    if (x$format==1){
        Pb206Pb204 <- x$x[,'Pb206Pb204']
        Pb207Pb204 <- x$x[,'Pb207Pb204']
        errPb206Pb204 <- x$x[,'errPb206Pb204']
        errPb207Pb204 <- x$x[,'errPb207Pb204']
        rho <- x$x[,'rhoXY']
    } else if (x$format==2){
        Pb206Pb204 <- 1/x$x[,'Pb204Pb206']
        Pb207Pb204 <- x$x[,'Pb207Pb206']/x$x[,'Pb204Pb206']
        J11 <- -1/x$x[,'Pb204Pb206']^2
        J22 <-  1/x$x[,'Pb204Pb206']
        J12 <- 0
        J21 <- -x$x[,'Pb207Pb206']/x$x[,'Pb204Pb206']^2
        E11 <- x$x[,'errPb204Pb206']^2
        E22 <- x$x[,'errPb207Pb206']^2
        E12 <- x$x[,'rhoXY']*x$x[,'errPb204Pb206']*x$x[,'errPb207Pb206']
        covmat <- errorprop(J11,J12,J21,J22,E11,E12,E22)
        errPb206Pb204 <- sqrt(covmat[,'varX'])
        errPb207Pb204 <- sqrt(covmat[,'varY'])
        rho <- covmat[,'cov']/(errPb206Pb204*errPb207Pb204)
    } else if (x$format==3){
        Pb206Pb204 <- 1/x$x[,'Pb204Pb206']
        Pb207Pb204 <- 1/x$x[,'Pb204Pb207']
        J11 <- -1/x$x[,'Pb204Pb206']^2
        J22 <- -1/x$x[,'Pb204Pb207']^2
        J12 <- 0
        J21 <- 0
        E11 <- x$x[,'errPb204Pb206']^2
        E22 <- x$x[,'errPb204Pb207']^2
        E12 <- x$x[,'rhoXY']*x$x[,'errPb204Pb206']*x$x[,'errPb204Pb207']
        covmat <- errorprop(J11,J12,J21,J22,E11,E12,E22)
        errPb206Pb204 <- sqrt(covmat[,'varX'])
        errPb207Pb204 <- sqrt(covmat[,'varY'])
        rho <- get.cor.div(Pb207Pb204,errPb207Pb204,
                           Pb206Pb204,errPb206Pb204,                            
                           Pb207Pb206,errPb207Pb206)
    }
    out <- cbind(Pb206Pb204,errPb206Pb204,Pb207Pb204,errPb207Pb204,rho)
    colnames(out) <- c('Pb206Pb204','errPb206Pb204',
                       'Pb207Pb204','errPb207Pb204','rho')
    out
}
# returns ratios, errors and correlations
PbPb.inverse.ratios <- function(x){
    ns <- length(x)
    out <- matrix(0,ns,5)
    if (x$format<2){
        Pb204Pb206 <- 1/x$x[,'Pb206Pb204']
        Pb207Pb206 <- x$x[,'Pb207Pb204']/x$x[,'Pb206Pb204']
        J11 <- -1/x$x[,'Pb206Pb204']^2
        J22 <-  1/x$x[,'Pb206Pb204']
        J12 <-  0
        J21 <- -x$x[,'Pb207Pb204']/x$x[,'Pb206Pb204']^2
        E11 <- x$x[,'Pb206Pb204']^2
        E22 <- x$x[,'Pb207Pb204']^2
        E12 <- x$x[,'rhoXY']*x$x[,'errPb206Pb204']*x$x[,'errPb207Pb204']
        covmat <- errorprop(J11,J12,J21,J22,E11,E12,E22)
        errPb204Pb206 <- sqrt(covmat[,'varX'])
        errPb207Pb206 <- sqrt(covmat[,'varY'])
        rho <- covmat[,'cov']/(errPb204Pb206*errPb207Pb206)
    } else {
        Pb204Pb206 <- x$x[,'Pb204Pb206']
        errPb204Pb206 <- x$x[,'errPb204Pb206']
        Pb207Pb206 <- x$x[,'Pb207Pb206']
        errPb207Pb206 <- x$x[,'errPb207Pb206']
        rho <- x$x[,'rhoXY']
    }
    out <- cbind(Pb204Pb206,errPb204Pb206,Pb207Pb206,errPb207Pb206,rho)
    colnames(out) <- c('Pb204Pb206','errPb204Pb206',
                       'Pb207Pb206','errPb207Pb206','rho')
    out
}

get.PbPb.age <- function(Pb207Pb206,sPb207Pb206,exterr=TRUE){
    c(0,0)
}
