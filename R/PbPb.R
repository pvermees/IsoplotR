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
        covmat <- errorprop(J11,J12,J21,J22,E11,E22,E12)
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
        covmat <- errorprop(J11,J12,J21,J22,E11,E22,E12)
        errPb206Pb204 <- x$x[,'errPb206Pb204']
        errPb207Pb204 <- x$x[,'errPb207Pb204']
        rho <- covmat[,'cov']/(errPb206Pb204*errPb207Pb204)
    }
    out <- cbind(Pb206Pb204,errPb206Pb204,Pb207Pb204,errPb207Pb204,rho)
    colnames(out) <- c('Pb206Pb204','errPb206Pb204',
                       'Pb207Pb204','errPb207Pb204','rho')
    out[abs(out[,'rho'])>1,'rho'] <- 0.99 # protect against rounding errors
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
        covmat <- errorprop(J11,J12,J21,J22,E11,E22,E12)
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
        covmat <- errorprop(J11,J12,J21,J22,E11,E22,E12)
        errPb204Pb206 <- sqrt(covmat[,'varX'])
        errPb207Pb206 <- x$x[,'errPb207Pb206']
        rho <- covmat[,'cov']/(errPb204Pb206*errPb207Pb206)
    }
    out <- cbind(Pb204Pb206,errPb204Pb206,Pb207Pb206,errPb207Pb206,rho)
    colnames(out) <- c('Pb204Pb206','errPb204Pb206',
                       'Pb207Pb206','errPb207Pb206','rho')
    out[abs(out[,'rho'])>1,'rho'] <- 0.99 # protect against rounding errors
    out
}

PbPb.age <- function(x,exterr=TRUE,i=NA,sigdig=NA,common.Pb=0){
    if (common.Pb == 0){
        y <- data2york(x,inverse=FALSE)
        PbPb <- quotient(y[,'X'],y[,'sX'],y[,'Y'],y[,'sY'],y[,'rXY'])
    } else if (common.Pb == 1){
        y <- data2york(x,inverse=FALSE)
        tt <- rep(1000,length(x))
        for (j in 1:10){
            i6474 <- stacey.kramers(tt)
            r64 <- y[,'X'] - i6474[,1]
            r74 <- y[,'Y'] - i6474[,2]
            PbPb <- quotient(r64,y[,'sX'],r74,y[,'sY'],y[,'rXY'])
            tt <- PbPb2t(PbPb)
        }
    } else if (common.Pb == 2){
        y <- data2york(x,inverse=TRUE)
        fit <- regression(y,model=1)
        PbPb <- get.76(y,a=fit$a[1],b=fit$b[1])
    } else if (common.Pb == 3){
        y <- data2york(x,inverse=FALSE)
        r64 <- y[,'X'] - settings('iratio','Pb206Pb204')[1]
        r74 <- y[,'Y'] - settings('iratio','Pb207Pb204')[1]
        PbPb <- quotient(r64,y[,'sX'],r74,y[,'sY'],y[,'rXY'])
        # alternative implementation:
        # y <- data2york(x,inverse=TRUE)
        # i74 <- settings('iratio','Pb207Pb204')[1]
        # i64 <- settings('iratio','Pb206Pb204')[1]
        # PbPb <- get.76(y,a=i74/i64,b=i74)
    }
    PbPb2t(PbPb,exterr=exterr,sigdig=sigdig,i=i)
}

PbPb2t <- function(PbPb,exterr=FALSE,sigdig=NA,i=NA){
    ns <- nrow(PbPb)
    out <- matrix(0,ns,2)
    colnames(out) <- c('t','s[t]')
    for (j in 1:ns){
        tt <- get.Pb207Pb206.age(PbPb[j,1],PbPb[j,2],exterr=exterr)
        out[j,] <- roundit(tt[1],tt[2],sigdig=sigdig)
    }
    if (!is.na(i)) out <- out[i,]
    out
}

get.76 <- function(y,a,b){
    ns <- nrow(y)
    out <- matrix(0,ns,2)
    out[,1] <- 2*a + b*y[,'X'] - y[,'Y']
    J1 <- b
    J2 <- rep(-1,ns)
    E11 <- y[,'sX']^2
    E22 <- y[,'sY']^2
    E12 <- y[,'sX']*y[,'sY']*y[,'rXY']
    out[,2] <- errorprop1x2(J1,J2,E11,E22,E12)
    out
}
