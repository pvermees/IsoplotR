# returns ratios, errors and correlations
PbPb.normal.ratios <- function(x){
    ns <- length(x)
    out <- matrix(0,ns,5)
    if (x$format==1){
        Pb206Pb204 <- x$x[,'Pb206Pb204']
        Pb207Pb204 <- x$x[,'Pb207Pb204']
        errPb206Pb204 <- x$x[,'errPb206Pb204']
        errPb207Pb204 <- x$x[,'errPb207Pb204']
        rho <- x$x[,'rho']
    } else if (x$format==2){
        Pb206Pb204 <- 1/x$x[,'Pb204Pb206']
        Pb207Pb204 <- x$x[,'Pb207Pb206']/x$x[,'Pb204Pb206']
        J11 <- -1/x$x[,'Pb204Pb206']^2
        J22 <-  1/x$x[,'Pb204Pb206']
        J12 <- rep(0,ns)
        J21 <- -x$x[,'Pb207Pb206']/x$x[,'Pb204Pb206']^2
        E11 <- x$x[,'errPb204Pb206']^2
        E22 <- x$x[,'errPb207Pb206']^2
        E12 <- x$x[,'rho']*x$x[,'errPb204Pb206']*x$x[,'errPb207Pb206']
        covmat <- errorprop(J11,J12,J21,J22,E11,E12,E22)
        errPb206Pb204 <- sqrt(covmat[,'varX'])
        errPb207Pb204 <- sqrt(covmat[,'varY'])
        rho <- covmat[,'cov']/(errPb206Pb204*errPb207Pb204)
    } else if (x$format==3){
        Pb206Pb204 <- x$x[,'Pb206Pb204']
        Pb207Pb204 <- x$x[,'Pb207Pb204']
        J11 <- rep(1,ns)
        J22 <- rep(1,ns)
        J12 <- rep(0,ns)
        J21 <- rep(0,ns)
        E11 <- x$x[,'errPb206Pb204']^2
        E22 <- x$x[,'errPb207Pb204']^2
        E12 <- get.cov.div(x$x[,'Pb207Pb204'],x$x[,'errPb207Pb204'],
                           x$x[,'Pb206Pb204'],x$x[,'errPb206Pb204'],
                           x$x[,'Pb207Pb206'],x$x[,'errPb207Pb206'])
        covmat <- errorprop(J11,J12,J21,J22,E11,E12,E22)
        errPb206Pb204 <- x$x[,'errPb206Pb204']
        errPb207Pb204 <- x$x[,'errPb207Pb204']
        rho <- covmat[,'cov']/(errPb206Pb204*errPb207Pb204)
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
    if (x$format==1){
        Pb204Pb206 <- 1/x$x[,'Pb206Pb204']
        Pb207Pb206 <- x$x[,'Pb207Pb204']/x$x[,'Pb206Pb204']
        J11 <- -1/x$x[,'Pb206Pb204']^2
        J22 <-  1/x$x[,'Pb206Pb204']
        J12 <-  rep(0,ns)
        J21 <- -x$x[,'Pb207Pb204']/x$x[,'Pb206Pb204']^2
        E11 <- x$x[,'errPb206Pb204']^2
        E22 <- x$x[,'errPb207Pb204']^2
        E12 <- x$x[,'rho']*x$x[,'errPb206Pb204']*x$x[,'errPb207Pb204']
        covmat <- errorprop(J11,J12,J21,J22,E11,E12,E22)
        errPb204Pb206 <- sqrt(covmat[,'varX'])
        errPb207Pb206 <- sqrt(covmat[,'varY'])
        rho <- covmat[,'cov']/(errPb204Pb206*errPb207Pb206)
    } else if (x$format==2){
        Pb204Pb206 <- x$x[,'Pb204Pb206']
        errPb204Pb206 <- x$x[,'errPb204Pb206']
        Pb207Pb206 <- x$x[,'Pb207Pb206']
        errPb207Pb206 <- x$x[,'errPb207Pb206']
        rho <- x$x[,'rho']
    } else if (x$format==3){
        Pb204Pb206 <- 1/x$x[,'Pb206Pb204']
        Pb207Pb206 <- x$x[,'Pb207Pb206']
        J11 <- -1/x$x[,'Pb206Pb204']^2
        J22 <- rep(1,ns)
        J12 <- rep(0,ns)
        J21 <- rep(0,ns)
        E11 <- x$x[,'errPb206Pb204']^2
        E22 <- x$x[,'errPb207Pb206']^2
        E12 <- get.cov.mult(x$x[,'Pb206Pb204'],x$x[,'errPb206Pb204'],
                            x$x[,'Pb207Pb206'],x$x[,'errPb207Pb206'],
                            x$x[,'Pb207Pb204'],x$x[,'errPb207Pb204'])
        covmat <- errorprop(J11,J12,J21,J22,E11,E12,E22)
        errPb204Pb206 <- sqrt(covmat[,'varX'])
        errPb207Pb206 <- x$x[,'errPb207Pb206']
        rho <- covmat[,'cov']/(errPb204Pb206*errPb207Pb206)
    }
    out <- cbind(Pb204Pb206,errPb204Pb206,Pb207Pb206,errPb207Pb206,rho)
    colnames(out) <- c('Pb204Pb206','errPb204Pb206',
                       'Pb207Pb206','errPb207Pb206','rho')
    out
}

PbPb.age <- function(x,exterr=TRUE,i=NA,sigdig=NA,i2i=TRUE,...){
    ns <- length(x)
    dat <- data2york(x,inverse=FALSE)
    if (i2i){
        fit <- isochron(x,plot=FALSE,exterr=exterr)        
        dat[,'Y'] <- dat[,'Y'] - fit$a[1]
        if (exterr) dat[,'sY'] <- sqrt(dat[,'sY']^2 + fit$a[2]^2)
    }
    out <- matrix(0,ns,2)
    colnames(out) <- c('t','s[t]')
    E <- matrix(0,2,2)
    J <- matrix(0,1,2)
    for (j in 1:ns) {
        E[1,1] <- dat[j,'sX']^2
        E[2,2] <- dat[j,'sY']^2
        E[1,2] <- dat[j,'rXY']*dat[j,'sX']*dat[j,'sY']
        E[2,1] <- E[1,2]
        J[1,1] <- -dat[j,'Y']/dat[j,'X']^2
        J[1,2] <- 1/dat[j,'X']        
        DP <- dat[j,'Y']/dat[j,'X']
        sDP <- sqrt(J%*%E%*%t(J))
        tt <- get.Pb207Pb206.age(DP,sDP,exterr=exterr)
        out[j,] <- roundit(tt[1],tt[2],sigdig=sigdig)
    }
    if (!is.na(i)) out <- out[i,]
    out
}
