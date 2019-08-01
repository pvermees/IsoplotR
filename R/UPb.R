wetherill <- function(x,i){
    out <- list()
    if (x$format < 4) labels <- c('Pb207U235','Pb206U238')
    else if (x$format < 7) labels <- c('Pb207U235','Pb206U238','Pb204U238')
    else labels <- c('Pb207U235','Pb206U238','Pb208Th232','Th232U238')
    if (x$format %in% c(1,3)){
        out$x <- x$x[i,labels]
        out$cov <- cor2cov2(x$x[i,'errPb207U235'],x$x[i,'errPb206U238'],x$x[i,'rhoXY'])
    } else if (x$format == 2){
        U238U235 <- iratio('U238U235')[1]
        Pb207U235 <- U238U235*x$x[i,'Pb207Pb206']/x$x[i,'U238Pb206']
        Pb206U238 <- 1/x$x[i,'U238Pb206']
        out$x <- c(Pb207U235,Pb206U238)
        J <- matrix(0,2,2)
        J[1,1] <- -Pb207U235/x$x[i,'U238Pb206']
        J[1,2] <- U238U235/x$x[i,'U238Pb206']
        J[2,1] <- -1/x$x[i,'U238Pb206']^2
        E <- cor2cov2(x$x[i,'errU238Pb206'],x$x[i,'errPb207Pb206'],x$x[i,'rhoXY'])
        out$cov <- J %*% E %*% t(J)
    } else if (x$format == 4){
        out$x <- x$x[i,labels]
        out$cov <- cor2cov3(x$x[i,'errPb207U235'],x$x[i,'errPb206U238'],
                            x$x[i,'errPb204U238'],x$x[i,'rhoXY'],
                            x$x[i,'rhoXZ'],x$x[i,'rhoYZ'])
    } else if (x$format == 5){
        U238U235 <- iratio('U238U235')[1]
        Pb207U235 <- U238U235*x$x[i,'Pb207Pb206']/x$x[i,'U238Pb206']
        Pb206U238 <- 1/x$x[i,'U238Pb206']
        Pb204U238 <- x$x[i,'Pb204Pb206']/x$x[i,'U238Pb206']
        out$x <- c(Pb207U235,Pb206U238,Pb204U238)
        J <- matrix(0,3,3)
        J[1,1] <- -Pb207U235/x$x[i,'U238Pb206']
        J[1,2] <- U238U235/x$x[i,'U238Pb206']
        J[2,1] <- -Pb206U238/x$x[i,'U238Pb206']
        J[3,1] <- -Pb204U238/x$x[i,'U238Pb206']
        J[3,3] <- 1/x$x[i,'U238Pb206']
        E <- cor2cov3(x$x[i,'errU238Pb206'],x$x[i,'errPb207Pb206'],
                      x$x[i,'errPb204Pb206'],x$x[i,'rhoXY'],
                      x$x[i,'rhoXZ'],x$x[i,'rhoYZ'])
        out$cov <- J %*% E %*% t(J)        
    } else if (x$format == 6){
        out$x <- x$x[i,labels]
        out$cov <- matrix(0,3,3)
        diag(out$cov) <-
            x$x[i,c('errPb207U235','errPb206U238','errPb204U238')]^2
        out$cov[1,2] <-
            get.cov.75.68(x$x[i,'Pb207U235'],x$x[i,'errPb207U235'],
                          x$x[i,'Pb206U238'],x$x[i,'errPb206U238'],
                          x$x[i,'Pb207Pb206'],x$x[i,'errPb207Pb206'])
        out$cov[1,3] <-
            get.cov.75.48(x$x[i,'Pb207U235'],x$x[i,'errPb207U235'],
                          x$x[i,'Pb204U238'],x$x[i,'errPb204U238'],
                          x$x[i,'Pb204Pb207'],x$x[i,'errPb204Pb207'])
        out$cov[2,3] <-
            get.cov.68.48(x$x[i,'Pb206U238'],x$x[i,'errPb206U238'],
                          x$x[i,'Pb204U238'],x$x[i,'errPb204U238'],
                          x$x[i,'Pb204Pb206'],x$x[i,'errPb204Pb206'])
        out$cov[2,1] <- out$cov[1,2]
        out$cov[3,1] <- out$cov[1,3]
        out$cov[3,2] <- out$cov[2,3]
    } else if (x$format == 7){
        out$x <- x$x[i,labels]
        out$cov <- cor2cov4(x$x[i,'errPb207U235'],x$x[i,'errPb206U238'],
                            x$x[i,'errPb208Th232'],x$x[i,'errTh232U238'],
                            x$x[i,'rhoXY'],x$x[i,'rhoXZ'],x$x[i,'rhoXW'],
                            x$x[i,'rhoYZ'],x$x[i,'rhoYW'],x$x[i,'rhoZW'])
    } else if (x$format == 8){
        U238U235 <- iratio('U238U235')[1]
        Pb207U235 <- U238U235*x$x[i,'Pb207Pb206']/x$x[i,'U238Pb206']
        Pb206U238 <- 1/x$x[i,'U238Pb206']
        Pb208Th232 <- x$x[i,'Pb208Pb206']/(x$x[i,'U238Pb206']*x$x[i,'Th232U238'])
        Th232U238 <- x$x[i,'Th232U238']
        out$x <- c(Pb207U235,Pb206U238,Pb208Th232,Th232U238)
        J <- matrix(0,4,4)
        J[1,1] <- -Pb207U235/x$x[i,'U238Pb206']
        J[1,2] <- U238U235/x$x[i,'U238Pb206']
        J[2,1] <- -Pb206U238/x$x[i,'U238Pb206']
        J[3,1] <- -Pb208Th232/x$x[i,'U238Pb206']
        J[3,3] <- 1/(x$x[i,'U238Pb206']*x$x[i,'Th232U238'])
        J[3,4] <- -Pb208Th232/x$x[i,'Th232U238']
        J[4,4] <- 1
        E <- cor2cov4(x$x[i,'errU238Pb206'],x$x[i,'errPb207Pb206'],
                      x$x[i,'errPb208Pb206'],x$x[i,'errTh232U238'],
                      x$x[i,'rhoXY'],x$x[i,'rhoXZ'],x$x[i,'rhoXW'],
                      x$x[i,'rhoYZ'],x$x[i,'rhoYW'],x$x[i,'rhoZW'])
        out$cov <- J %*% E %*% t(J)        
    }
    names(out$x) <- labels
    colnames(out$cov) <- labels
    rownames(out$cov) <- labels
    class(out) <- "wetherill"
    out
}
tera.wasserburg <- function(x,i){
    out <- list()
    if (x$format < 4) labels <- c('U238Pb206','Pb207Pb206')
    else if (x$format < 7) labels <- c('U238Pb206','Pb207Pb206','Pb204Pb206')
    else labels <- c('U238Pb206','Pb207Pb206','Pb208Pb206','Th232U238')
    if (x$format==1){
        U238U235 <- iratio('U238U235')[1]
        U238Pb206 <- 1/x$x[i,'Pb206U238']
        Pb207Pb206 <- x$x[i,'Pb207U235']/(x$x[i,'Pb206U238']*U238U235)
        J <- matrix(0,2,2)
        J[1,2] <- -1/x$x[i,'Pb206U238']^2
        J[2,1] <- 1/(x$x[i,'Pb206U238']*U238U235)
        J[2,2] <- -Pb207Pb206/x$x[i,'Pb206U238']
        E <- cor2cov2(x$x[i,'errPb207U235'],x$x[i,'errPb206U238'],x$x[i,'rhoXY'])
        out$x <- c(U238Pb206,Pb207Pb206)
        out$cov <- J %*% E %*% t(J)
    } else if (x$format == 2){
        out$x <- x$x[i,labels]
        out$cov <- matrix(0,2,2)
        diag(out$cov) <- x$x[i,c('errU238Pb206','errPb207Pb206')]^2
        out$cov <- cor2cov2(x$x[i,'errU238Pb206'],x$x[i,'errPb207Pb206'],x$x[i,'rhoXY'])
    } else if (x$format == 3){
        U238Pb206 <- 1/x$x[i,'Pb206U238']
        out$x <- c(U238Pb206,x$x[i,'Pb207Pb206'])
        J <- matrix(0,2,2)
        E <- cor2cov2(x$x[i,'errPb206U238'],x$x[i,'errPb207Pb206'],x$x[i,'rhoYZ'])
        J[1,1] <- -U238Pb206^2
        J[2,2] <- 1
        out$cov <- J %*% E %*% t(J)
    } else if (x$format == 4){
        U238U235 <- iratio('U238U235')[1]
        U238Pb206 <- 1/x$x[i,'Pb206U238']
        Pb207Pb206 <- x$x[i,'Pb207U235']/(x$x[i,'Pb206U238']*U238U235)
        Pb204Pb206 <- x$x[i,'Pb204U238']/x$x[i,'Pb206U238']
        E <- cor2cov3(x$x[i,'errPb207U235'],x$x[i,'errPb206U238'],
                      x$x[i,'errPb204U238'],x$x[i,'rhoXY'],
                      x$x[i,'rhoXZ'],x$x[i,'rhoYZ'])
        J <- matrix(0,3,3)
        J[1,2] <- -U238Pb206/x$x[i,'Pb206U238']
        J[2,1] <- 1/(x$x[i,'Pb206U238']*U238U235)
        J[2,2] <- -Pb207Pb206/x$x[i,'Pb206U238']
        J[3,2] <- -Pb204Pb206/x$x[i,'Pb206U238']
        J[3,3] <- 1/x$x[i,'Pb206U238']
        out$x <- c(U238Pb206,Pb207Pb206,Pb204Pb206)
        out$cov <- J %*% E %*% t(J)
    } else if (x$format == 5){
        out$x <- x$x[i,labels]
        out$cov <- matrix(0,3,3)
        diag(out$cov) <- x$x[i,c('errU238Pb206','errPb207Pb206','errPb204Pb206')]^2
        out$cov <-
            cor2cov3(x$x[i,'errU238Pb206'],x$x[i,'errPb207Pb206'],x$x[i,'errPb204Pb206'],
                     x$x[i,'rhoXY'],x$x[i,'rhoXZ'],x$x[i,'rhoYZ'])
    } else if (x$format == 6){
        U238Pb206 <- 1/x$x[i,'Pb206U238']
        Pb207Pb206 <- x$x[i,'Pb207Pb206']
        Pb204Pb206 <- x$x[i,'Pb204Pb206']
        out$x <- c(U238Pb206,Pb207Pb206,Pb204Pb206)
        out$cov <- matrix(0,3,3)
        out$cov[1,1] <- (U238Pb206*x$x[i,'errPb206U238']/x$x[i,'Pb206U238'])^2
        out$cov[2,2] <- x$x[i,'errPb207Pb206']^2
        out$cov[3,3] <- x$x[i,'errPb204Pb206']^2
        out$cov[1,2] <- get.cov.76.86(x$x[i,'Pb207Pb206'],x$x[i,'errPb207Pb206'],
                                      x$x[i,'Pb206U238'],x$x[i,'errPb206U238'],
                                      x$x[i,'Pb207U235'],x$x[i,'errPb207U235'])
        out$cov[1,3] <- get.cov.46.86(x$x[i,'Pb204Pb206'],x$x[i,'errPb204Pb206'],
                                      x$x[i,'Pb206U238'],x$x[i,'errPb206U238'],
                                      x$x[i,'Pb204U238'],x$x[i,'errPb204U238'])
        out$cov[2,3] <- get.cov.46.76(x$x[i,'Pb204Pb206'],x$x[i,'errPb204Pb206'],
                                      x$x[i,'Pb207Pb206'],x$x[i,'errPb207Pb206'],
                                      x$x[i,'Pb204Pb207'],x$x[i,'errPb204Pb207'])
        out$cov[2,1] <- out$cov[1,2]
        out$cov[3,1] <- out$cov[1,3]
        out$cov[3,2] <- out$cov[2,3]
    } else if (x$format == 7){
        U238U235 <- iratio('U238U235')[1]
        U238Pb206 <- 1/x$x[i,'Pb206U238']
        Pb207Pb206 <- x$x[i,'Pb207U235']/(x$x[i,'Pb206U238']*U238U235)
        Pb208Pb206 <- x$x[i,'Pb208Th232']*x$x[i,'Th232U238']/x$x[i,'Pb206U238']
        Th232U238 <- x$x[i,'Th232U238']
        E <- cor2cov4(x$x[i,'errPb207U235'],x$x[i,'errPb206U238'],
                      x$x[i,'errPb208Th232'],x$x[i,'errTh232U238'],
                      x$x[i,'rhoXY'],x$x[i,'rhoXZ'],x$x[i,'rhoXW'],
                      x$x[i,'rhoYZ'],x$x[i,'rhoYW'],x$x[i,'rhoZW'])
        J <- matrix(0,4,4)
        J[1,2] <- -U238Pb206/x$x[i,'Pb206U238']
        J[2,1] <- 1/(x$x[i,'Pb206U238']*U238U235)
        J[2,2] <- -Pb207Pb206/x$x[i,'Pb206U238']
        J[3,2] <- -Pb208Pb206/x$x[i,'Pb206U238']
        J[3,3] <- x$x[i,'Th232U238']/x$x[i,'Pb206U238']
        J[3,4] <- x$x[i,'Pb208Th232']/x$x[i,'Pb206U238']
        J[4,4] <- 1
        out$x <- c(U238Pb206,Pb207Pb206,Pb208Pb206,Th232U238)
        out$cov <- J %*% E %*% t(J)
    } else if (x$format == 8){
        out$x <- x$x[i,labels]
        out$cov <- cor2cov4(x$x[i,'errU238Pb206'],x$x[i,'errPb207Pb206'],
                            x$x[i,'errPb208Pb206'],x$x[i,'errTh232U238'],
                            x$x[i,'rhoXY'],x$x[i,'rhoXZ'],x$x[i,'rhoXW'],
                            x$x[i,'rhoYZ'],x$x[i,'rhoYW'],x$x[i,'rhoZW'])
    }
    names(out$x) <- labels
    colnames(out$cov) <- labels
    rownames(out$cov) <- labels
    class(out) <- "terawasserburg"
    out
}
get.UPb.isochron.ratios.204 <- function(x,i){
    if (x$format%in%c(4,5,6)){
        labels <- c('U238Pb206','Pb204Pb206','U235Pb207','Pb204Pb207')
    } else {
        stop('Incorrect input format for the get.UPb.isochron.ratios function.')
    }
    U <- iratio('U238U235')[1]
    tw <- tera.wasserburg(x,i) # 38/06, 07/06 and 04/06
    U8Pb6 <- tw$x['U238Pb206']
    Pb46 <- tw$x['Pb204Pb206']
    U5Pb7 <- tw$x['U238Pb206']/(U*tw$x['Pb207Pb206'])
    Pb47 <- tw$x['Pb204Pb206']/tw$x['Pb207Pb206']
    J <- matrix(0,4,3)
    J[1,1] <- 1
    J[2,2] <- 1
    J[3,1] <- 1/(U*tw$x['Pb207Pb206'])
    J[3,2] <- -U5Pb7/tw$x['Pb207Pb206']
    J[4,2] <- -Pb47/tw$x['Pb207Pb206']
    J[4,3] <- 1/tw$x['Pb207Pb206']
    out <- list()
    out$x <- c(U8Pb6,Pb46,U5Pb7,Pb47)
    out$cov <- J %*% tw$cov %*% t(J)
    names(out$x) <- labels
    colnames(out$cov) <- labels
    rownames(out$cov) <- labels
    out
}
get.UPb.isochron.ratios.208 <- function(x,i,tt=0){
    if (x$format%in%c(7,8)){
        labels <- c('U238Pb206','Pb208cPb206','U235Pb207',
                    'Pb208cPb207','Th232U238','Th232Pb208')
    } else {
        stop('Incorrect input format for the get.UPb.isochron.ratios function.')
    }
    l2 <- settings('lambda','Th232')[1]
    U <- iratio('U238U235')[1]
    tw <- tera.wasserburg(x,i) # 38/06, 07/06, 08/06, 32/38
    U8Pb6 <- tw$x['U238Pb206']
    Pb8c6 <- tw$x['Pb208Pb206'] - tw$x['Th232U238']*tw$x['U238Pb206']*(exp(l2*tt)-1)
    U5Pb7 <- tw$x['U238Pb206']/(U*tw$x['Pb207Pb206'])
    Pb8c7 <- Pb8c6/tw$x['Pb207Pb206']
    Th2Pb8 <- tw$x['Th232U238']*tw$x['U238Pb206']/tw$x['Pb208Pb206']
    J <- matrix(0,6,4)
    J[1,1] <- 1
    J[2,1] <- -tw$x['Th232U238']*(exp(l2*tt)-1)
    J[2,3] <- 1
    J[2,4] <- -tw$x['U238Pb206']*(exp(l2*tt)-1)
    J[3,1] <- 1/(U*tw$x['Pb207Pb206'])
    J[3,2] <- -U5Pb7/tw$x['Pb207Pb206']
    J[4,1] <- J[2,1]/tw$x['Pb207Pb206']
    J[4,2] <- -Pb8c7/tw$x['Pb207Pb206']
    J[4,3] <- J[2,3]/tw$x['Pb207Pb206']
    J[4,4] <- J[2,4]/tw$x['Pb207Pb206']
    J[5,4] <- 1
    J[6,1] <- tw$x['Th232U238']/tw$x['Pb208Pb206']
    J[6,3] <- -Th2Pb8/tw$x['Pb208Pb206']
    J[6,4] <- tw$x['U238Pb206']/tw$x['Pb208Pb206']
    out <- list()
    out$x <- c(U8Pb6,Pb8c6,U5Pb7,Pb8c7,tw$x['Th232U238'],Th2Pb8)
    out$cov <- J %*% tw$cov %*% t(J)
    names(out$x) <- labels
    colnames(out$cov) <- labels
    rownames(out$cov) <- labels
    out
}

w2tw <- function(w){
    U <- iratio('U238U235')[1]
    Pb207U235 <- w[,1]
    errPb207U235 <- w[,2]
    Pb206U238 <- w[,3]
    errPb206U238 <- w[,4]
    rho <- w[,5]
    U238Pb206 <- 1/Pb206U238
    Pb207Pb206 <- Pb207U235/(U*Pb206U238)
    J11 <- 0*U238Pb206
    J12 <- -U238Pb206/Pb206U238
    J21 <- 1/(U*Pb206U238)
    J22 <- -Pb207Pb206/Pb206U238
    E11 <- errPb207U235^2
    E22 <- errPb206U238^2
    E12 <- rho*errPb207U235*errPb206U238
    err <- errorprop(J11,J12,J21,J22,E11,E22,E12)
    errU238Pb206 <- sqrt(err[,'varX'])
    errPb207Pb206 <- sqrt(err[,'varY'])
    rho <- err[,'cov']/(errU238Pb206*errPb207Pb206)
    out <- cbind(U238Pb206,errU238Pb206,Pb207Pb206,errPb207Pb206,rho)
    colnames(out) <- c('U238Pb206','errU238Pb206',
                       'Pb207Pb206','errPb207Pb206','rhoXY')
    out
}
tw2w <- function(tw){
    U <- iratio('U238U235')[1]
    U238Pb206 <- tw[,1]
    errU238Pb206 <- tw[,2]
    Pb207Pb206 <- tw[,3]
    errPb207Pb206 <- tw[,4]
    rho <- tw[,5]
    Pb207U235 <- U*Pb207Pb206/U238Pb206
    Pb206U238 <- 1/U238Pb206
    J11 <- -Pb207U235/U238Pb206
    J12 <- U/U238Pb206
    J21 <- -1/U238Pb206^2
    J22 <- 0*Pb206U238
    E11 <- errU238Pb206^2
    E22 <- errPb207Pb206^2
    E12 <- rho*errU238Pb206*errPb207Pb206
    err <- errorprop(J11,J12,J21,J22,E11,E22,E12)
    errPb207U235 <- sqrt(err[,'varX'])
    errPb206U238 <- sqrt(err[,'varY'])
    rho <- err[,'cov']/(errPb207U235*errPb206U238)
    out <- cbind(Pb207U235,errPb207U235,Pb206U238,errPb206U238,rho)
    colnames(out) <- c('Pb207U235','errPb207U235',
                       'Pb206U238','errPb206U238','rhoXY')
    out
}

# convert data to a 5-column table for concordia analysis
flat.UPb.table <- function(x,type=1){
    ns <- length(x)
    out <- matrix(0,ns,5)
    if (type==1){
        colnames(out) <- c('Pb207U235','errPb207U235',
                           'Pb206U238','errPb206U238','rhoXY')
    } else if (type==2){
        colnames(out) <- c('U238Pb206','errU238Pb206',
                           'Pb207Pb206','errPb207Pb206','rhoXY')
    } else if (type==3){
        colnames(out) <- c('Pb206U238','errPb206U238',
                           'Pb208Th232','errPb208Th232','rhoXY')
    } else {
        stop('Incorrect concordia type.')
    }
    for (i in 1:ns){
        if (type==1){
            xi <- wetherill(x,i)
            j1 <- 'Pb207U235'
            j2 <- 'Pb206U238'
        } else if (type==2){
            xi <- tera.wasserburg(x,i)
            j1 <- 'U238Pb206'
            j2 <- 'Pb207Pb206'
        } else if (type==3){
            xi <- wetherill(x,i)
            j1 <- 'Pb206U238'
            j2 <- 'Pb208Th232'
        }
        out[i,1] <- xi$x[j1]
        out[i,3] <- xi$x[j2]
        out[i,2] <- sqrt(xi$cov[j1,j1])
        out[i,4] <- sqrt(xi$cov[j2,j2])
        out[i,5] <- xi$cov[j1,j2]/(out[i,2]*out[i,4])
    }
    out
}

age_to_wetherill_ratios <- function(tt,st=0,exterr=FALSE,d=diseq()){
    out <- list()
    labels <- c('Pb207U235','Pb206U238')
    l8 <- settings('lambda','U238')[1]
    l5 <- settings('lambda','U235')[1]
    D <- mclean(tt=tt,d=d,exterr=exterr)
    out$x <- c(D$Pb207U235,D$Pb206U238)
    E <- matrix(0,3,3)
    diag(E) <- c(st,lambda('U235')[2],lambda('U238')[2])^2
    J <- matrix(0,2,3)
    J[1,1] <- D$dPb207U235dt
    J[2,1] <- D$dPb206U238dt
    J[1,2] <- D$dPb207U235dl35
    J[2,3] <- D$dPb206U238dl38
    out$cov <- J %*% E %*% t(J)
    names(out$x) <- labels
    rownames(out$cov) <- labels
    colnames(out$cov) <- labels
    out
}
age_to_terawasserburg_ratios <- function(tt,st=0,exterr=FALSE,d=diseq()){
    out <- list()
    labels <- c('U238Pb206','Pb207Pb206')
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    U <- iratio('U238U235')[1]
    tt <- check.zero.UPb(tt)
    D <- mclean(tt=tt,d=d,exterr=exterr)
    U238Pb206 <- 1/D$Pb206U238
    Pb207Pb206 <- D$Pb207U235/(U*D$Pb206U238)
    d75dt <- D$dPb207U235dt
    d68dt <- D$dPb206U238dt
    d75dl35 <- D$dPb207U235dl35
    d68dl38 <- D$dPb206U238dl38
    out$x <- c(U238Pb206,Pb207Pb206)
    E <- matrix(0,4,4)
    diag(E) <- c(st,lambda('U235')[2],lambda('U238')[2],iratio('U238U235')[2])^2
    J <- matrix(0,2,4)
    J[1,1] <- -d68dt/D$Pb206U238^2
    J[2,1] <- (d75dt*D$Pb206U238-D$Pb207U235*d68dt)/(U*D$Pb206U238^2)
    J[1,3] <- -d68dl38/D$Pb206U238^2
    J[2,2] <- d75dl35/(U*D$Pb206U238)
    J[2,3] <- -Pb207Pb206*d68dl38/D$Pb206U238
    J[2,4] <- -Pb207Pb206/U
    out$cov <- J %*% E %*% t(J)
    names(out$x) <- labels
    rownames(out$cov) <- labels
    colnames(out$cov) <- labels
    out
}
age_to_cottle_ratios <- function(tt,st=0,exterr=FALSE,d=diseq()){
    out <- list()
    labels <- c('Pb206U238','Pb208Th232')
    l8 <- settings('lambda','U238')[1]
    l2 <- settings('lambda','Th232')[1]
    D <- mclean(tt=tt,d=d,exterr=exterr)
    Pb6U8 <- D$Pb206U238
    Pb8Th2 <- exp(l2*tt)-1
    out$x <- c(Pb6U8,Pb8Th2)
    E <- matrix(0,3,3)
    diag(E) <- c(st,lambda('U238')[2],lambda('Th232')[2])^2
    J <- matrix(0,2,3)
    J[1,1] <- D$dPb206U238dt
    J[2,1] <- l2*exp(l2*tt)
    J[1,2] <- D$dPb206U238dl38
    J[2,3] <- tt*exp(l2*tt)
    out$cov <- J %*% E %*% t(J)
    names(out$x) <- labels
    rownames(out$cov) <- labels
    colnames(out$cov) <- labels
    out
}

age_to_Pb207U235_ratio <- function(tt,st=0,d=diseq()){
    ns <- length(tt)
    if (ns>1){
        out <- matrix(0,ns,2)
        colnames(out) <- c('75','s[75]')
        for (i in 1:ns){
            if (length(st)<ns) sti <- st[1]
            else sti <- st[i]
            out[i,] <- age_to_Pb207U235_ratio(tt[i],st=sti,d=d[i])
        }
    } else {
        l5 <- lambda('U235')[1]
        D <- mclean(tt=tt,d=d)
        R <- D$Pb207U235
        J <- D$dPb207U235dt
        R.err <- abs(J*st)
        out <- cbind(R,R.err)
        colnames(out) <- c('75','s[75]')
    }
    out
}
age_to_U235Pb207_ratio <- function(tt,st=0,d=diseq()){
    ns <- length(tt)
    if (ns>1){
        out <- matrix(0,ns,2)
        colnames(out) <- c('57','s[57]')
        for (i in 1:ns){
            if (length(st)<ns) sti <- st[1]
            else sti <- st[i]
            out[i,] <- age_to_U235Pb207_ratio(tt[i],st=sti,d=d[i])
        }
    } else {
        l5 <- lambda('U235')[1]
        D <- mclean(tt=tt,d=d)
        R <- 1/D$Pb207U235
        J <- -R*D$dPb207U235dt/D$Pb207U235
        R.err <- abs(J*st)
        out <- cbind(R,R.err)
        colnames(out) <- c('57','s[57]')
    }
    out
}
age_to_Pb206U238_ratio <- function(tt,st=0,d=diseq()){
    ns <- length(tt)
    if (ns>1){
        out <- matrix(0,ns,2)
        colnames(out) <- c('68','s[68]')
        for (i in 1:ns){
            if (length(st)<ns) sti <- st[1]
            else sti <- st[i]
            out[i,] <- age_to_Pb206U238_ratio(tt[i],st=sti,d=d[i])
        }
    } else {
        l8 <- lambda('U238')[1]
        D <- mclean(tt=tt,d=d)
        R <- D$Pb206U238
        J <- D$dPb206U238dt
        R.err <- abs(J*st)
        out <- cbind(R,R.err)
        colnames(out) <- c('68','s[68]')
    }
    out
}
age_to_U238Pb206_ratio <- function(tt,st=0,d=diseq()){
    ns <- length(tt)
    if (ns>1){
        out <- matrix(0,ns,2)
        colnames(out) <- c('86','s[86]')
        for (i in 1:ns){
            if (length(st)<ns) sti <- st[1]
            else sti <- st[i]
            out[i,] <- age_to_U238Pb206_ratio(tt[i],st=sti,d=d[i])
        }
    } else {
        l8 <- lambda('U238')[1]
        tt <- check.zero.UPb(tt)
        D <- mclean(tt=tt,d=d)
        R <- 1/D$Pb206U238
        J <- -D$dPb206U238dt/D$Pb206U238^2
        R.err <- abs(J*st)
        out <- cbind(R,R.err)
        colnames(out) <- c('86','s[86]')
    }
    out
}
age_to_Pb207Pb206_ratio <- function(tt,st=0,d=diseq()){
    ns <- length(tt)
    if (ns>1){
        out <- matrix(0,ns,2)
        colnames(out) <- c('76','s[76]')
        for (i in 1:ns){
            if (length(st)<ns) sti <- st[1]
            else sti <- st[i]
            out[i,] <- age_to_Pb207Pb206_ratio(tt[i],st=sti,d=d[i])
        }
    } else {
        l5 <- lambda('U235')[1]
        l8 <- lambda('U238')[1]
        tt <- check.zero.UPb(tt)
        U <- iratio('U238U235')[1]
        D <- mclean(tt=tt,d=d)
        R <- (1/U)*D$Pb207U235/D$Pb206U238
        J <- (1/U)*(D$dPb207U235dt*D$Pb206U238 -
                    D$Pb207U235*D$dPb206U238dt)/D$Pb206U238^2
        R.err <- abs(J*st)
        out <- cbind(R,R.err)
        colnames(out) <- c('76','s[76]')
    }
    out
}
age_to_Pb208Th232_ratio <- function(tt,st=0){
    ns <- length(tt)
    if (ns>1){
        out <- matrix(0,ns,2)
        colnames(out) <- c('82','s[82]')
        for (i in 1:ns){
            if (length(st)<ns) sti <- st[1]
            else sti <- st[i]
            out[i,] <- age_to_Pb208Th232_ratio(tt[i],st=sti)
        }
    } else {
        l2 <- lambda('Th232')[1]
        R <- exp(l2*tt)-1
        J <- l2*exp(l2*tt)
        R.err <- abs(J*st)
        out <- cbind(R,R.err)
        colnames(out) <- c('82','s[82]')
    }
    out    
}

check.zero.UPb <- function(tt){
    smallnum <- 2*.Machine$double.neg.eps/lambda('U238')[1]
    if (length(tt)>1){
        pos <- which(tt>0)
        out <- tt
        out[!pos] <- smallnum
    } else {
        out <- max(tt,smallnum)
    }
    out
}

get.Pb204U238.ratios <- function(x){
    labels <- c('Pb204U238','errPb204U238')
    if (x$format %in% c(1,2,3,7,8)){
        stop('No 204Pb measurements available!')
    } else if (x$format %in% c(4,6)){
        out <- subset(x$x,select=labels)
    } else {
        Pb204U238 <- x$x[,'Pb204Pb206']/x$x[,'U238Pb206']
        errPb204U238 <-
            Pb204U238*sqrt( (x$x[,'errPb204Pb206']/x$x[,'Pb204Pb206'])^2 +
                            (x$x[,'errU238Pb206']/x$x[,'U238Pb206'])^2 )
        out <- cbind(Pb204U238,errPb204U238)
    }
    colnames(out) <- labels
    out
}
get.Pb207U235.ratios <- function(x,exterr=FALSE){
    ns <- length(x)
    out <- matrix(0,ns,2)
    labels <- c('Pb207U235','errPb207U235')
    colnames(out) <- labels
    if (x$format %in% c(1,3,4,6,7)){
        out <- subset(x$x,select=labels)
    } else if (x$format %in% c(2,5,8)){
        R <- iratio('U238U235')[1]
        sR <- iratio('U238U235')[2]
        X <- x$x[,'U238Pb206']
        sX <- x$x[,'errU238Pb206']
        Y <- x$x[,'Pb207Pb206']
        sY <- x$x[,'errPb207Pb206']
        covXY <- x$x[,'rhoXY']*sX*sY
        out[,'Pb207U235'] <- R*Y/X
        relerr2 <- (sX/X)^2 -2*covXY/(X*Y) + (sY/Y)^2
        if (exterr) relerr2 <- relerr2 + (sR/R)^2
        out[,'errPb207U235'] <- sqrt(relerr2)*out[,'Pb207U235']
    }
    out
}
get.U235Pb207.ratios <- function(x,exterr=FALSE){
    ns <- length(x)
    out <- matrix(0,ns,2)
    labels <- c('U235Pb207','errU235Pb207')
    colnames(out) <- labels
    if (x$format %in% c(1,3,4,6,7)){
        out[,'U235Pb207'] <- 1/x$x[,'Pb207U235']
        out[,'errU235Pb207'] <- out[,'U235Pb207']*
            x$x[,'errPb207U235']/x$x[,'Pb207U235']
    } else if (x$format %in% c(2,5,8)){
        R <- iratio('U238U235')[1]
        sR <- iratio('U238U235')[2]
        X <- x$x[,'U238Pb206']
        sX <- x$x[,'errU238Pb206']
        Y <- x$x[,'Pb207Pb206']
        sY <- x$x[,'errPb207Pb206']
        covXY <- x$x[,'rhoXY']*sX*sY
        out[,'U235Pb207'] <- X/(R*Y)
        relerr2 <- (sX/X)^2 -2*covXY/(X*Y) + (sY/Y)^2
        if (exterr) relerr2 <- relerr2 + (sR/R)^2
        out[,'errU235Pb207'] <- sqrt(relerr2)*out[,'U235Pb207']
    }
    out
}
get.Pb206U238.ratios <- function(x){
    ns <- length(x)
    out <- matrix(0,ns,2)
    labels <- c('Pb206U238','errPb206U238')
    colnames(out) <- labels
    if (x$format %in% c(1,3,4,6,7)){
        out <- subset(x$x,select=labels)
    } else if (x$format %in% c(2,5,8)){
        out[,'Pb206U238'] <- 1/x$x[,'U238Pb206']
        out[,'errPb206U238'] <- out[,'Pb206U238']*
            x$x[,'errU238Pb206']/x$x[,'U238Pb206']
    }
    out
}
get.U238Pb206.ratios <- function(x){
    ns <- length(x)
    out <- matrix(0,ns,2)
    labels <- c('U238Pb206','errU238Pb206')
    colnames(out) <- labels
    if (x$format %in% c(1,3,4,6,7)){
        out[,'U238Pb206'] <- 1/x$x[,'Pb206U238']
        out[,'errU238Pb206'] <- out[,'U238Pb206']*
            x$x[,'errPb206U238']/x$x[,'Pb206U238']
    } else if (x$format %in% c(2,5,8)){
        out <- subset(x$x,select=labels)
    }
    out
}
get.Pb207Pb206.ratios <- function(x,exterr=FALSE){
    ns <- length(x)
    out <- matrix(0,ns,2)
    labels <- c('Pb207Pb206','errPb207Pb206')
    colnames(out) <- labels
    if (x$format %in% c(1,4,7)){
        R <- iratio('U238U235')[1]
        sR <- iratio('U238U235')[2]
        X <- x$x[,'Pb207U235']
        sX <- x$x[,'errPb207U235']
        Y <- x$x[,'Pb206U238']
        sY <- x$x[,'errPb206U238']
        rhoXY <- x$x[,'rhoXY']
        covXY <- rhoXY*sX*sY
        out[,'Pb207Pb206'] <- X/(R*Y)
        relerr2 <- (sX/X)^2 - 2*covXY/(X*Y) + (sY/Y)^2
        if (exterr) relerr2 <- relerr2 + (sR/R)^2
        out[,'errPb207Pb206'] <- sqrt(relerr2)*out[,'Pb207Pb206']
    } else if (x$format %in% c(2,3,5,6,8)){
        out <- subset(x$x,select=labels)
    }
    out
}
get.Pb208Th232.ratios <- function(x){
    ns <- length(x)
    out <- matrix(0,ns,2)
    labels <- c('Pb208Th232','errPb208Th232')
    colnames(out) <- labels
    if (x$format == 7){
        out <- subset(x$x,select=labels)
    } else if (x$format == 8){
        out[,'Pb208Th232'] <- x$x[,'Pb208Pb206']/(x$x[,'U238Pb206']*x$x[,'Th232U238'])
        J1 <- -out[,'Pb208Th232']/x$x[,'U238Pb206']
        J2 <- 1/(x$x[,'U238Pb206']*x$x[,'Th232U238'])
        J3 <- -out[,'Pb208Th232']/x$x[,'Th232U238']
        E11 <- x$x[,'errU238Pb206']^2
        E22 <- x$x[,'errPb208Pb206']^2
        E33 <- x$x[,'errTh232U238']^2
        E12 <- x$x[,'rhoXZ']*x$x[,'errU238Pb206']*x$x[,'errPb208Pb206']
        E13 <- x$x[,'rhoXW']*x$x[,'errU238Pb206']*x$x[,'errTh232U238']
        E23 <- x$x[,'rhoZW']*x$x[,'errPb208Pb206']*x$x[,'errTh232U238']
        out[,'errPb208Th232'] <- errorprop1x3(J1,J2,J3,E11,E22,E33,E12,E13,E23)
    } else {
        stop('Wrong input format: no Pb208 or Th232 present in this dataset.')
    }
    out
}
get.Pb208Pb206.ratios <- function(x){
    ns <- length(x)
    out <- matrix(0,ns,2)
    labels <- c('Pb208Pb206','errPb208Pb206')
    colnames(out) <- labels
    if (x$format == 7){
        out[,'Pb208Pb206'] <- x$x[,'Pb208Th232']*x$x[,'Th232U238']/x$x[,'Pb206U238']
        J1 <- -out[,'Pb208Pb206']/x$x[,'Pb206U238'] # d/dPb6U8 = d/dY
        J2 <- x$x[,'Th232U238']/x$x[,'Pb206U238']   # d/dPb8Th2 = d/dZ
        J3 <- x$x[,'Pb208Th232']/x$x[,'Pb206U238']  # d/dTh2U8 = d/dW
        sX <- x$x[,'errPb206U238']
        sY <- x$x[,'errPb208Th232']
        sZ <- x$x[,'errTh232U238']
        rXY <- x$x[,'rhoYZ']
        rXZ <- x$x[,'rhoYW']
        rYZ <- x$x[,'rhoZW']
        E12 <- rXY*sX*sY
        E13 <- rXZ*sX*sZ
        E23 <- rYZ*sY*sZ
        out[,'errPb208Pb206'] <- errorprop1x3(J1,J2,J3,sX^2,sY^2,sZ^2,E12,E13,E23)
    } else if (x$format == 8){
        out <- x$x[,labels]
    } else {
        stop('Wrong input format: no Pb208 present in this dataset.')
    }
    out
}

get.Pb207U235.age <- function(x,...){ UseMethod("get.Pb207U235.age",x) }
get.Pb207U235.age.default <- function(x,sx=0,exterr=FALSE,d=diseq(),...){
    ns <- length(x)
    if (ns>1){
        out <- matrix(0,ns,2)
        colnames(out) <- c('t75','s[t75]')
        for (i in 1:ns){
            if (length(sx) < ns) sxi <- sx[1]
            else sxi <- sx[i]
            out[i,] <- get.Pb207U235.age(x[i],sxi,exterr=exterr,d=d[i])
        }
    } else {
        l5 <- lambda('U235')[1]
        sl5 <- lambda('U235')[2]
        if (x>-1) t.75 <- log(1+x)/l5 else t.75 <- 0
        J <- matrix(0,1,2)
        if (d$equilibrium){
            J[1,1] <- 1/(l5*(1+x))                       # dt/dx
            if (exterr & x>-1) J[1,2] <- log(1+x)/l5^2   # dt/dl5
        } else { # apply a disequilibrium correction
            search.range <- c(t.75/1000,t.75+100)
            t.75 <- stats::optimize(diseq.75.misfit,interval=search.range,x=x,d=d)$minimum
            D <- mclean(tt=t.75,d=d,exterr=exterr)    # implicit differentiation of 
            J[1,1] <- -1/D$dPb207U235dt               # mf=(x-Pb7U5)^2 => dt/dx
            J[1,2] <- D$dPb207U235dl35/D$dPb207U235dt # and dt/dl35
        }
        E <- matrix(0,2,2)
        E[1,1] <- sx^2
        E[2,2] <- sl5^2
        st.75 <- sqrt(J %*% E %*% t(J))
        out <- c(t.75,st.75)
        names(out) <- c('t75','s[t75]')
    }
    out
}
get.Pb207U235.age.UPb <- function(x,i,exterr=FALSE,...){
    r75 <- get.Pb207U235.ratios(x)
    get.Pb207U235.age(r75[i,'Pb207U235'],r75[i,'errPb207U235'],exterr=exterr,d=x$d)
}
get.Pb207U235.age.wetherill <- function(x,exterr=FALSE,...){
    i <- 'Pb207U235'
    r75 <- x$x[i]
    sr75 <- sqrt(x$cov[i,i])
    get.Pb207U235.age(r75,sr75,exterr=exterr,d=x$d,...)
}

get.Pb206U238.age <- function(x,...){ UseMethod("get.Pb206U238.age",x) }
get.Pb206U238.age.default <- function(x,sx=0,exterr=FALSE,d=diseq(),...){
    ns <- length(x)
    if (ns>1){
        out <- matrix(0,ns,2)
        colnames(out) <- c('t68','s[t68]')
        for (i in 1:ns){
            if (length(sx) < ns) sxi <- sx[1]
            else sxi <- sx[i]
            out[i,] <- get.Pb206U238.age(x[i],sxi,exterr=exterr,d=d[i])
        }
    } else {
        l8 <- lambda('U238')[1]
        sl8 <- lambda('U238')[2]
        if (x>-1) t.init <- log(1+x)/l8 else t.init <- 0
        J <- matrix(0,1,2)
        if (d$equilibrium){
            t.68 <- t.init
            J[1,1] <- 1/(l8*(1+x))                       # dt/dx
            if (exterr & x>-1) J[1,2] <- log(1+x)/l8^2   # dt/dl38
        } else { # apply a disequilibrium correction
            t.68 <- tryCatch({
                search.range <- c(t.init/1000,t.init+100)
                stats::optimise(diseq.68.misfit,interval=search.range,x=x,d=d)$minimum
            }, error = function(error_condition) {
                t.init
            })
            D <- mclean(tt=t.68,d=d,exterr=exterr)    # implicit differentiation of 
            J[1,1] <- -1/D$dPb206U238dt               # mf=(x-Pb6U8)^2 => dt/dx
            J[1,2] <- D$dPb206U238dl38/D$dPb206U238dt # and dt/dl38
        }
        E <- matrix(0,2,2)
        E[1,1] <- sx^2
        E[2,2] <- sl8^2
        st.68 <- sqrt(J %*% E %*% t(J))
        out <- c(t.68,st.68)
        names(out) <- c('t68','s[t68]')
    }
    out
}
get.Pb206U238.age.UPb <- function(x,i=NA,exterr=FALSE,...){
    r68 <- get.Pb206U238.ratios(x)
    if (is.na(i)){
        ns <- length(x)
        out <- matrix(0,ns,2)
        for (j in 1:ns){
            out[j,] <- get.Pb206U238.age.UPb(x,i=j,exterr=exterr,...)
        }
    } else {
        out <- get.Pb206U238.age(r68[i,'Pb206U238'],r68[i,'errPb206U238'],
                                 exterr=exterr,d=x$d[i],...)
    }
    out
}
get.Pb206U238.age.wetherill <- function(x,exterr=FALSE,...){
    i <- 'Pb206U238'
    r68 <- x$x[i]
    sr68 <- sqrt(x$cov[i,i])
    get.Pb206U238.age(r68,sr68,exterr=exterr,d=x$d,...)
}
get.Pb206U238.age.terawasserburg <- function(x,exterr=FALSE,...){
    i <- 'U238Pb206'
    r86 <- x$x[i]
    r68 <- 1/r86
    sr68 <- sqrt(x$cov[i,i])/r86
    get.Pb206U238.age(r68,sr68,exterr=exterr,d=x$d,...)
}

twslope <- function(tt=0,d=diseq()){
    D <- mclean(tt=tt,d=d)
    D$dPb207U235dt/D$dPb206U238dt
}

get.Pb207Pb206.age <- function(x,...){ UseMethod("get.Pb207Pb206.age",x) }
get.Pb207Pb206.age.default <- function(x,sx=0,exterr=FALSE,d=diseq(),t.68=NA,...){
    ns <- length(x)
    if (ns>1){
        out <- matrix(0,ns,2)
        colnames(out) <- c('t76','s[t76]')
        for (i in 1:ns){
            if (length(sx) < ns) sxi <- sx[1]
            else sxi <- sx[i]
            out[i,] <- get.Pb207Pb206.age(x[i],sxi,exterr=exterr,d=d[i],t.68=t.68)
        }
    } else {
        interval <- c(1/10000,10000)
        if (!d$equilibrium & !any(is.na(t.68))){
            midpoint <- stats::optimise(twslope,d=d,interval=interval)$minimum
            if (t.68<midpoint){
                interval[2] <- midpoint
            } else {
                interval[1] <- midpoint
            }
        }
        t.76 <- stats::optimise(get.76.misfit,x=x,d=d,interval=interval)$minimum
        D <- mclean(tt=t.76,d=d,exterr=exterr)
        J <- matrix(0,1,4)
        J[1,1] <- -1/D$dPb207Pb206dt                # dt/dx
        J[1,2] <- D$dPb207Pb206dl35/D$dPb207Pb206dt # dt/dl35
        J[1,3] <- D$dPb207Pb206dl38/D$dPb207Pb206dt # dt/dl38
        J[1,4] <- D$dPb207Pb206dU/D$dPb207Pb206dt   # dt/dU
        E <- matrix(0,4,4)
        E[1,1] <- sx^2
        E[2,2] <- lambda('U235')[2]^2
        E[3,3] <- lambda('U238')[2]^2
        E[4,4] <- iratio('U238U235')[2]^2
        st.76 <- sqrt( J %*% E %*% t(J) )
        out <- c(t.76,st.76)
        names(out) <- c('t76','s[t76]')
    }
    out    
}
get.Pb207Pb206.age.UPb <- function(x,i,exterr=FALSE,...){
    r76 <- get.Pb207Pb206.ratios(x)
    get.Pb207Pb206.age(r76[i,'Pb207Pb206'],r76[i,'errPb207Pb206'],
                       exterr=exterr,d=x$d[i],...)
}
get.Pb207Pb206.age.wetherill <- function(x,exterr=FALSE,...){
    U <- iratio('U238U235')[1]
    r76 <- x$x['Pb207U235']/(U*x$x['Pb206U238'])
    J <- matrix(0,1,2)
    E <- x$cov
    J[1,1] <- 1/(U*x$x['Pb206U238'])                   # d76d75
    J[1,2] <- -x$x['Pb207U235']/(U*x$x['Pb206U238']^2) # d76d68
    sr76 <- J %*% E %*% t(J)
    get.Pb207Pb206.age(r76,sr76,exterr=exterr,d=x$d)
}
get.Pb207Pb206.age.terawasserburg <- function(x,exterr=FALSE,...){
    r76 <- x$x['Pb207Pb206']
    sr76 <- x$cov['Pb207Pb206','Pb207Pb206']
    get.Pb207Pb206.age(r76,sr76,exterr=exterr,d=x$d)
}

get.Pb208Th232.age <- function(x,...){ UseMethod("get.Pb208Th232.age",x) }
get.Pb208Th232.age.default <- function(x,sx=0,exterr=FALSE,...){
    ns <- length(x)
    if (ns>1){
        out <- matrix(0,ns,2)
        colnames(out) <- c('t82','s[t82]')
        for (i in 1:ns){
            if (length(sx) < ns) sxi <- sx[1]
            else sxi <- sx[i]
            out[i,] <- get.Pb208Th232.age(x[i],sxi,exterr=exterr)
        }
    } else {
        l2 <- lambda('Th232')[1]
        sl2 <- lambda('Th232')[2]
        if (x>-1) t.init <- log(1+x)/l2 else t.init <- 0
        J <- matrix(0,1,2)
        t.82 <- t.init
        J[1,1] <- 1/(l2*(1+x))                       # dt/dx
        if (exterr & x>-1) J[1,2] <- log(1+x)/l2^2   # dt/dl2
        E <- matrix(0,2,2)
        E[1,1] <- sx^2
        E[2,2] <- sl2^2
        st.82 <- sqrt(J %*% E %*% t(J))
        out <- c(t.82,st.82)
        names(out) <- c('t82','s[t82]')
    }
    out
}
get.Pb208Th232.age.UPb <- function(x,i=NA,exterr=FALSE,...){
    r82 <- get.Pb208Th232.ratios(x)
    if (is.na(i)){
        ns <- length(x)
        out <- matrix(0,ns,2)
        for (j in 1:ns){
            out[j,] <- get.Pb208Th232.age.UPb(x,i=j,exterr=exterr,...)
        }
    } else {
        out <- get.Pb208Th232.age(r82[i,'Pb208Th232'],r82[i,'errPb208Th232'],
                                 exterr=exterr,...)
    }
    out
}

# x is an object of class \code{UPb}
# returns a matrix of 7/5, 6/8, 7/6
# and concordia ages and their uncertainties.
UPb.age <- function(x,exterr=FALSE,i=NA,sigdig=NA,conc=TRUE,show.p=FALSE,common.Pb=0,...){
    if (common.Pb>0) X <- Pb0corr(x,option=common.Pb)
    else X <- x
    labels <- c('t.75','s[t.75]','t.68','s[t.68]','t.76','s[t.76]')
    hasTh <- x$format%in%c(7,8)
    if (hasTh) labels <- c(labels,'t.82','s[t.82]')
    if (conc) labels <- c(labels,'t.conc','s[t.conc]')
    if (conc & show.p) labels <- c(labels,'p[conc]')
    if (!is.na(i)){
        t.75 <- get.Pb207U235.age(X,i,exterr=exterr)
        t.68 <- get.Pb206U238.age(X,i,exterr=exterr)
        t.76 <- get.Pb207Pb206.age(X,i,exterr=exterr,t.68=t.68[1])
        t.75.out <- roundit(t.75[1],t.75[2],sigdig=sigdig)
        t.68.out <- roundit(t.68[1],t.68[2],sigdig=sigdig)
        t.76.out <- roundit(t.76[1],t.76[2],sigdig=sigdig)
        out <- c(t.75.out,t.68.out,t.76.out)
        if (hasTh){
            t.82 <- get.Pb208Th232.age(X,i,exterr=exterr)
            t.82.out <- roundit(t.82[1],t.82[2],sigdig=sigdig)
            out <- c(out,t.82.out)
        }
        if (conc){
            t.conc <- concordia.age(x=X,i=i,exterr=exterr)
            t.conc.out <- roundit(t.conc$age[1],t.conc$age[2],sigdig=sigdig)
            out <- c(out,t.conc.out)
        }
        if (conc & show.p){
            SS.concordance <-
                LL.concordia.age(tt=t.conc$age[1],cc=wetherill(X,i),
                                 mswd=TRUE,exterr=exterr,d=x$d)
            p.value <- 1-stats::pchisq(SS.concordance,1)
            if (!is.na(sigdig)) p.value <- signif(p.value,sigdig)
            out <- c(out,p.value)
        }
        names(out) <- labels
    } else {
        nn <- nrow(X$x)
        out <- matrix(0,nn,length(labels))
        for (i in 1:nn){
            out[i,] <- UPb.age(X,i=i,exterr=exterr,sigdig=sigdig,
                               conc=conc,show.p=show.p)
        }
        colnames(out) <- labels
    }
    out
}

filter.UPb.ages <- function(x,type=4,cutoff.76=1100,exterr=FALSE,
                            cutoff.disc=list(-15,5,TRUE),common.Pb=0){
    tout <- UPb.age(x,exterr=exterr,conc=(type==5),common.Pb=common.Pb)
    if (any(is.na(cutoff.disc))){
        is.concordant <- rep(TRUE,length(x))
    } else {
        # apply cutoff filter before common Pb correction?
        if (cutoff.disc[[3]] & common.Pb>0){
            tin <- UPb.age(x,exterr=exterr,conc=(type==5),common.Pb=0)
        } else { # apply cutoff filter after common Pb correction.
            tin <- tout
        }
        is.concordant <- concordant(tin,cutoff.76=cutoff.76,cutoff.disc=cutoff.disc)
        if (!any(is.concordant)){
            stop(paste0('There are no concordant grains in this sample.',
                        'Try adjusting the discordance limits OR ',
                        'apply a common-Pb correction OR ',
                        '(if you have already applied a common-Pb correction), ',
                        'apply the discordance filter before the ',
                        'common-Pb correction.'))
        }
    }
    out <- matrix(NA,length(x),2)
    if (type==1){
        out[is.concordant,] <- tout[is.concordant,c('t.75','s[t.75]'),drop=FALSE]
    } else if (type==2){
        out[is.concordant,] <- tout[is.concordant,c('t.68','s[t.68]'),drop=FALSE]
    } else if (type==3){
        out[is.concordant,] <- tout[is.concordant,c('t.76','s[t.76]'),drop=FALSE]
    } else if (type==4){
        do.76 <- (tout[,'t.68']>cutoff.76)
        i.76 <- as.vector(which(do.76 & is.concordant))
        i.68 <- as.vector(which(!do.76 & is.concordant))
        out[i.76,] <- tout[i.76,c('t.76','s[t.76]'),drop=FALSE]
        out[i.68,] <- tout[i.68,c('t.68','s[t.68]'),drop=FALSE]
    } else if (type==5){
        out[is.concordant,] <- tout[is.concordant,c('t.conc','s[t.conc]'),drop=FALSE]
    } else if (type==6){
        out[is.concordant,] <- tout[is.concordant,c('t.82','s[t.82]'),drop=FALSE]
    }
    colnames(out) <- c('t','s[t]')
    out
}

concordant <- function(tt,cutoff.76=1100,cutoff.disc=list(-15,5,TRUE)){
    do.76 <- (tt[,'t.68'] > cutoff.76)
    disc.75.68 <- 100*(1-tt[,'t.75']/tt[,'t.68'])
    disc.68.76 <- 100*(1-tt[,'t.68']/tt[,'t.76'])
    conc.75.68 <- !do.76 & (disc.75.68>cutoff.disc[[1]]) & (disc.75.68<cutoff.disc[[2]])
    conc.68.76 <- do.76 & (disc.68.76>cutoff.disc[[1]]) & (disc.68.76<cutoff.disc[[2]])
    conc.75.68 | conc.68.76
}
