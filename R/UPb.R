wetherill <- function(x,i=1,format){
    if (missing(format)){
        X <- x$x
        format <- x$format
    } else {
        X <- x
    }
    out <- list()
    if (format < 4) labels <- c('Pb207U235','Pb206U238')
    else if (format < 7) labels <- c('Pb207U235','Pb206U238','Pb204U238')
    else labels <- c('Pb207U235','Pb206U238','Pb208Th232','Th232U238')
    if (format %in% c(1,3)){
        out$x <- X[i,labels]
        out$cov <- cor2cov2(X[i,'errPb207U235'],X[i,'errPb206U238'],X[i,'rhoXY'])
    } else if (format == 2){
        U238U235 <- iratio('U238U235')[1]
        Pb207U235 <- U238U235*X[i,'Pb207Pb206']/X[i,'U238Pb206']
        Pb206U238 <- 1/X[i,'U238Pb206']
        out$x <- c(Pb207U235,Pb206U238)
        J <- matrix(0,2,2)
        J[1,1] <- -Pb207U235/X[i,'U238Pb206']
        J[1,2] <- U238U235/X[i,'U238Pb206']
        J[2,1] <- -1/X[i,'U238Pb206']^2
        E <- cor2cov2(X[i,'errU238Pb206'],X[i,'errPb207Pb206'],X[i,'rhoXY'])
        out$cov <- J %*% E %*% t(J)
    } else if (format == 4){
        out$x <- X[i,labels]
        out$cov <- cor2cov3(X[i,'errPb207U235'],X[i,'errPb206U238'],
                            X[i,'errPb204U238'],X[i,'rhoXY'],
                            X[i,'rhoXZ'],X[i,'rhoYZ'])
    } else if (format == 5){
        U238U235 <- iratio('U238U235')[1]
        Pb207U235 <- U238U235*X[i,'Pb207Pb206']/X[i,'U238Pb206']
        Pb206U238 <- 1/X[i,'U238Pb206']
        Pb204U238 <- X[i,'Pb204Pb206']/X[i,'U238Pb206']
        out$x <- c(Pb207U235,Pb206U238,Pb204U238)
        J <- matrix(0,3,3)
        J[1,1] <- -Pb207U235/X[i,'U238Pb206']
        J[1,2] <- U238U235/X[i,'U238Pb206']
        J[2,1] <- -Pb206U238/X[i,'U238Pb206']
        J[3,1] <- -Pb204U238/X[i,'U238Pb206']
        J[3,3] <- 1/X[i,'U238Pb206']
        E <- cor2cov3(X[i,'errU238Pb206'],X[i,'errPb207Pb206'],
                      X[i,'errPb204Pb206'],X[i,'rhoXY'],
                      X[i,'rhoXZ'],X[i,'rhoYZ'])
        out$cov <- J %*% E %*% t(J)        
    } else if (format == 6){
        out$x <- X[i,labels]
        out$cov <- matrix(0,3,3)
        diag(out$cov) <-
            X[i,c('errPb207U235','errPb206U238','errPb204U238')]^2
        out$cov[1,2] <-
            get.cov.75.68(X[i,'Pb207U235'],X[i,'errPb207U235'],
                          X[i,'Pb206U238'],X[i,'errPb206U238'],
                          X[i,'Pb207Pb206'],X[i,'errPb207Pb206'])
        out$cov[1,3] <-
            get.cov.75.48(X[i,'Pb207U235'],X[i,'errPb207U235'],
                          X[i,'Pb204U238'],X[i,'errPb204U238'],
                          X[i,'Pb204Pb207'],X[i,'errPb204Pb207'])
        out$cov[2,3] <-
            get.cov.68.48(X[i,'Pb206U238'],X[i,'errPb206U238'],
                          X[i,'Pb204U238'],X[i,'errPb204U238'],
                          X[i,'Pb204Pb206'],X[i,'errPb204Pb206'])
        out$cov[2,1] <- out$cov[1,2]
        out$cov[3,1] <- out$cov[1,3]
        out$cov[3,2] <- out$cov[2,3]
    } else if (format == 7){
        out$x <- X[i,labels]
        out$cov <- cor2cov4(X[i,'errPb207U235'],X[i,'errPb206U238'],
                            X[i,'errPb208Th232'],X[i,'errTh232U238'],
                            X[i,'rhoXY'],X[i,'rhoXZ'],X[i,'rhoXW'],
                            X[i,'rhoYZ'],X[i,'rhoYW'],X[i,'rhoZW'])
    } else if (format == 8){
        U238U235 <- iratio('U238U235')[1]
        Pb207U235 <- U238U235*X[i,'Pb207Pb206']/X[i,'U238Pb206']
        Pb206U238 <- 1/X[i,'U238Pb206']
        Pb208Th232 <- X[i,'Pb208Pb206']/(X[i,'U238Pb206']*X[i,'Th232U238'])
        Th232U238 <- X[i,'Th232U238']
        out$x <- c(Pb207U235,Pb206U238,Pb208Th232,Th232U238)
        J <- matrix(0,4,4)
        J[1,1] <- -Pb207U235/X[i,'U238Pb206']
        J[1,2] <- U238U235/X[i,'U238Pb206']
        J[2,1] <- -Pb206U238/X[i,'U238Pb206']
        J[3,1] <- -Pb208Th232/X[i,'U238Pb206']
        J[3,3] <- 1/(X[i,'U238Pb206']*X[i,'Th232U238'])
        J[3,4] <- -Pb208Th232/X[i,'Th232U238']
        J[4,4] <- 1
        E <- cor2cov4(X[i,'errU238Pb206'],X[i,'errPb207Pb206'],
                      X[i,'errPb208Pb206'],X[i,'errTh232U238'],
                      X[i,'rhoXY'],X[i,'rhoXZ'],X[i,'rhoXW'],
                      X[i,'rhoYZ'],X[i,'rhoYW'],X[i,'rhoZW'])
        out$cov <- J %*% E %*% t(J)        
    }
    names(out$x) <- labels
    colnames(out$cov) <- labels
    rownames(out$cov) <- labels
    class(out) <- "wetherill"
    out
}
tera.wasserburg <- function(x,i=1,format){
    if (missing(format)){
        X <- x$x
        format <- x$format
    } else {
        X <- x
    }
    out <- list()
    if (format < 4) labels <- c('U238Pb206','Pb207Pb206')
    else if (format < 7) labels <- c('U238Pb206','Pb207Pb206','Pb204Pb206')
    else labels <- c('U238Pb206','Pb207Pb206','Pb208Pb206','Th232U238')
    if (format==1){
        U238U235 <- iratio('U238U235')[1]
        U238Pb206 <- 1/X[i,'Pb206U238']
        Pb207Pb206 <- X[i,'Pb207U235']/(X[i,'Pb206U238']*U238U235)
        J <- matrix(0,2,2)
        J[1,2] <- -1/X[i,'Pb206U238']^2
        J[2,1] <- 1/(X[i,'Pb206U238']*U238U235)
        J[2,2] <- -Pb207Pb206/X[i,'Pb206U238']
        E <- cor2cov2(X[i,'errPb207U235'],X[i,'errPb206U238'],X[i,'rhoXY'])
        out$x <- c(U238Pb206,Pb207Pb206)
        out$cov <- J %*% E %*% t(J)
    } else if (format == 2){
        out$x <- X[i,labels]
        out$cov <- matrix(0,2,2)
        diag(out$cov) <- X[i,c('errU238Pb206','errPb207Pb206')]^2
        out$cov <- cor2cov2(X[i,'errU238Pb206'],X[i,'errPb207Pb206'],X[i,'rhoXY'])
    } else if (format == 3){
        U238Pb206 <- 1/X[i,'Pb206U238']
        out$x <- c(U238Pb206,X[i,'Pb207Pb206'])
        J <- matrix(0,2,2)
        E <- cor2cov2(X[i,'errPb206U238'],X[i,'errPb207Pb206'],X[i,'rhoYZ'])
        J[1,1] <- -U238Pb206^2
        J[2,2] <- 1
        out$cov <- J %*% E %*% t(J)
    } else if (format == 4){
        U238U235 <- iratio('U238U235')[1]
        U238Pb206 <- 1/X[i,'Pb206U238']
        Pb207Pb206 <- X[i,'Pb207U235']/(X[i,'Pb206U238']*U238U235)
        Pb204Pb206 <- X[i,'Pb204U238']/X[i,'Pb206U238']
        E <- cor2cov3(X[i,'errPb207U235'],X[i,'errPb206U238'],
                      X[i,'errPb204U238'],X[i,'rhoXY'],
                      X[i,'rhoXZ'],X[i,'rhoYZ'])
        J <- matrix(0,3,3)
        J[1,2] <- -U238Pb206/X[i,'Pb206U238']
        J[2,1] <- 1/(X[i,'Pb206U238']*U238U235)
        J[2,2] <- -Pb207Pb206/X[i,'Pb206U238']
        J[3,2] <- -Pb204Pb206/X[i,'Pb206U238']
        J[3,3] <- 1/X[i,'Pb206U238']
        out$x <- c(U238Pb206,Pb207Pb206,Pb204Pb206)
        out$cov <- J %*% E %*% t(J)
    } else if (format == 5){
        out$x <- X[i,labels]
        out$cov <- matrix(0,3,3)
        diag(out$cov) <- X[i,c('errU238Pb206','errPb207Pb206','errPb204Pb206')]^2
        out$cov <-
            cor2cov3(X[i,'errU238Pb206'],X[i,'errPb207Pb206'],X[i,'errPb204Pb206'],
                     X[i,'rhoXY'],X[i,'rhoXZ'],X[i,'rhoYZ'])
    } else if (format == 6){
        U238Pb206 <- 1/X[i,'Pb206U238']
        Pb207Pb206 <- X[i,'Pb207Pb206']
        Pb204Pb206 <- X[i,'Pb204Pb206']
        out$x <- c(U238Pb206,Pb207Pb206,Pb204Pb206)
        out$cov <- matrix(0,3,3)
        out$cov[1,1] <- (U238Pb206*X[i,'errPb206U238']/X[i,'Pb206U238'])^2
        out$cov[2,2] <- X[i,'errPb207Pb206']^2
        out$cov[3,3] <- X[i,'errPb204Pb206']^2
        out$cov[1,2] <- get.cov.76.86(X[i,'Pb207Pb206'],X[i,'errPb207Pb206'],
                                      X[i,'Pb206U238'],X[i,'errPb206U238'],
                                      X[i,'Pb207U235'],X[i,'errPb207U235'])
        out$cov[1,3] <- get.cov.46.86(X[i,'Pb204Pb206'],X[i,'errPb204Pb206'],
                                      X[i,'Pb206U238'],X[i,'errPb206U238'],
                                      X[i,'Pb204U238'],X[i,'errPb204U238'])
        out$cov[2,3] <- get.cov.46.76(X[i,'Pb204Pb206'],X[i,'errPb204Pb206'],
                                      X[i,'Pb207Pb206'],X[i,'errPb207Pb206'],
                                      X[i,'Pb204Pb207'],X[i,'errPb204Pb207'])
        out$cov[2,1] <- out$cov[1,2]
        out$cov[3,1] <- out$cov[1,3]
        out$cov[3,2] <- out$cov[2,3]
    } else if (format == 7){
        U238U235 <- iratio('U238U235')[1]
        U238Pb206 <- 1/X[i,'Pb206U238']
        Pb207Pb206 <- X[i,'Pb207U235']/(X[i,'Pb206U238']*U238U235)
        Pb208Pb206 <- X[i,'Pb208Th232']*X[i,'Th232U238']/X[i,'Pb206U238']
        Th232U238 <- X[i,'Th232U238']
        E <- cor2cov4(X[i,'errPb207U235'],X[i,'errPb206U238'],
                      X[i,'errPb208Th232'],X[i,'errTh232U238'],
                      X[i,'rhoXY'],X[i,'rhoXZ'],X[i,'rhoXW'],
                      X[i,'rhoYZ'],X[i,'rhoYW'],X[i,'rhoZW'])
        J <- matrix(0,4,4)
        J[1,2] <- -U238Pb206/X[i,'Pb206U238']
        J[2,1] <- 1/(X[i,'Pb206U238']*U238U235)
        J[2,2] <- -Pb207Pb206/X[i,'Pb206U238']
        J[3,2] <- -Pb208Pb206/X[i,'Pb206U238']
        J[3,3] <- X[i,'Th232U238']/X[i,'Pb206U238']
        J[3,4] <- X[i,'Pb208Th232']/X[i,'Pb206U238']
        J[4,4] <- 1
        out$x <- c(U238Pb206,Pb207Pb206,Pb208Pb206,Th232U238)
        out$cov <- J %*% E %*% t(J)
    } else if (format == 8){
        out$x <- X[i,labels]
        out$cov <- cor2cov4(X[i,'errU238Pb206'],X[i,'errPb207Pb206'],
                            X[i,'errPb208Pb206'],X[i,'errTh232U238'],
                            X[i,'rhoXY'],X[i,'rhoXZ'],X[i,'rhoXW'],
                            X[i,'rhoYZ'],X[i,'rhoYW'],X[i,'rhoZW'])
    }
    names(out$x) <- labels
    colnames(out$cov) <- labels
    rownames(out$cov) <- labels
    class(out) <- "terawasserburg"
    out
}
get.UPb.isochron.ratios.204 <- function(x,i=NULL){
    if (x$format%in%c(4,5,6)){
        labels <- c('U238Pb206','Pb204Pb206','U235Pb207','Pb204Pb207')
    } else {
        stop('Format does not contain 204Pb.')
    }
    if (is.null(i)){
        ns <- length(x)
        out <- matrix(0,ns,length(labels))
        for (j in 1:ns){
            out[j,] <- get.UPb.isochron.ratios.204(x,i=j)$x
        }
        colnames(out) <- labels
        return(out)
    }
    U <- iratio('U238U235')[1]
    tw <- tera.wasserburg(x,i) # 38/06, 07/06 and 04/06
    U8Pb6 <- tw$x['U238Pb206']
    Pb46 <- tw$x['Pb204Pb206']
    U5Pb7 <- tw$x['U238Pb206']/(U*tw$x['Pb207Pb206'])
    Pb47 <- tw$x['Pb204Pb206']/tw$x['Pb207Pb206']
    J <- matrix(0,4,3)
    J[1,1] <- 1
    J[2,3] <- 1
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
get.UPb.isochron.ratios.208 <- function(x,i=NULL,tt=0){
    if (x$format>6){
        labels <- c('U238Pb206','Pb208cPb206','U235Pb207',
                    'Pb208cPb207','Th232U238','Th232Pb208',
                    'Pb206cPb208','Pb207cPb208')
    } else {
        stop('Incorrect input format for the get.UPb.isochron.ratios.208 function.')
    }
    if (is.null(i)){
        ns <- length(x)
        out <- matrix(0,ns,length(labels))
        for (j in 1:ns){
            out[j,] <- get.UPb.isochron.ratios.208(x,i=j,tt=tt)$x
        }
        colnames(out) <- labels
        return(out)
    }
    D <- mclean(tt,d=x$d[i])
    l2 <- settings('lambda','Th232')[1]
    U <- iratio('U238U235')[1]
    tw <- tera.wasserburg(x,i) # 38/06, 07/06, 08/06, 32/38
    U8Pb6 <- tw$x['U238Pb206']
    Pb8c6 <- tw$x['Pb208Pb206'] -
        tw$x['Th232U238']*tw$x['U238Pb206']*(exp(l2*tt)-1)
    U5Pb7 <- tw$x['U238Pb206']/(U*tw$x['Pb207Pb206'])
    Pb8c7 <- Pb8c6/tw$x['Pb207Pb206']
    Th2Pb8 <- tw$x['Th232U238']*tw$x['U238Pb206']/tw$x['Pb208Pb206']
    Pb6c8 <- 1/tw$x['Pb208Pb206'] -
        D$Pb206U238*tw$x['U238Pb206']/tw$x['Pb208Pb206']
    Pb7c8 <- tw$x['Pb207Pb206']/tw$x['Pb208Pb206'] -
        D$Pb207U235*tw$x['U238Pb206']/(U*tw$x['Pb208Pb206'])
    J <- matrix(0,8,4)
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
    J[7,1] <- -D$Pb206U238/tw$x['Pb208Pb206']
    J[7,3] <- -Pb6c8/tw$x['Pb208Pb206']
    J[8,1] <- -D$Pb207U235/(U*tw$x['Pb208Pb206'])
    J[8,2] <- 1/tw$x['Pb208Pb206']
    J[8,3] <- -Pb7c8/tw$x['Pb208Pb206']
    out <- list()
    out$x <- c(U8Pb6,Pb8c6,U5Pb7,Pb8c7,
               tw$x['Th232U238'],Th2Pb8,Pb6c8,Pb7c8)
    out$cov <- J %*% tw$cov %*% t(J)
    names(out$x) <- labels
    colnames(out$cov) <- labels
    rownames(out$cov) <- labels
    out
}

w2tw <- function(w,format){
    if (format %in% c(1,2,3)){
        cnames <- c('U238Pb206','errU238Pb206',
                    'Pb207Pb206','errPb207Pb206','rhoXY')
    } else if (format%in%c(4,5,6)){
        cnames <- c('U238Pb206','errU238Pb206',
                    'Pb207Pb206','errPb207Pb206',
                    'Pb206Pb206','errPb204Pb206',
                    'rhoXY','rhoXZ','rhoYZ')
    } else if (format%in%c(7,8)){
        cnames <- c('U238Pb206','errU238Pb206',
                    'Pb207Pb206','errPb207Pb206',
                    'Pb208Pb206','errPb208Pb206',
                    'Th232U238','errTh232U238',
                    'rhoXY','rhoXZ','rhoXW',
                    'rhoYZ','rhoYW','rhoZW')
    } else {
        stop('Invalid input format.')
    }
    ns <- nrow(w)
    out <- w*0
    colnames(out) <- cnames
    for (i in 1:ns){
        tw <- tera.wasserburg(x=w,i=i,format=format)
        out[i,] <- wtw_helper(x=tw$x,covmat=tw$cov,cnames=cnames)
    }
    out
}

tw2w <- function(tw,format){
    if (format %in% c(1,2,3)){
        cnames <- c('Pb207U235','errPb207U235',
                    'Pb206U238','errPb206U238','rhoXY')
    } else if (format%in%c(4,5,6)){
        cnames <- c('Pb207U235','errPb207U235',
                    'Pb206U238','errPb206U238',
                    'Pb204U238','errPb204U238',
                    'rhoXY','rhoXZ','rhoYZ')
    } else if (format%in%c(7,8)){
        cnames <- c('Pb207U235','errPb207U235',
                    'Pb206U238','errPb206U238',
                    'Pb208Th232','errPb208Th232',
                    'Th232U238','errTh232U238',
                    'rhoXY','rhoXZ','rhoXW',
                    'rhoYZ','rhoYW','rhoZW')
    } else {
        stop('Invalid input format.')
    }
    ns <- nrow(tw)
    out <- tw*0
    colnames(out) <- cnames
    for (i in 1:ns){
        w <- wetherill(x=tw,i=i,format=format)
        out[i,] <- wtw_helper(x=w$x,covmat=w$cov,cnames=cnames)
    }
    out
}

wtw_helper <- function(x,covmat,cnames){
    nc <- length(cnames)
    out <- rep(0,nc)
    names(out) <- cnames
    err <- sqrt(diag(covmat))
    cormat <- 0*covmat
    pos <- which(diag(covmat)>0)
    cormat[pos,pos] <- stats::cov2cor(covmat[pos,pos])
    out[c(1,3)] <- x[1:2]
    out[c(2,4)] <- err[1:2]
    out['rhoXY'] <- cormat[1,2]
    if (nc>5){
        out[5] <- x[3]
        out[6] <- err[3]
        out['rhoXZ'] <- cormat[1,3]
        out['rhoYZ'] <- cormat[2,3]
    }
    if (nc>9){
        out[7] <- x[4]
        out[8] <- err[4]
        out['rhoXW'] <- cormat[1,4]
        out['rhoYW'] <- cormat[2,4]
        out['rhoZW'] <- cormat[3,4]
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

age_to_Pb207U235_ratio <- function(tt,st=0,d=diseq(),exterr=FALSE){
    ns <- length(tt)
    if (ns>1){
        out <- matrix(0,ns,2)
        colnames(out) <- c('75','s[75]')
        for (i in 1:ns){
            if (length(st)<ns) sti <- st[1]
            else sti <- st[i]
            out[i,] <- age_to_Pb207U235_ratio(tt[i],st=sti,d=d[i],exterr=exterr)
        }
    } else {
        l5 <- lambda('U235')[1]
        D <- mclean(tt=tt,d=d,exterr=exterr)
        R <- D$Pb207U235
        J <- D$dPb207U235dt
        if (exterr){
            J <- c(J,D$dPb207U235dl35)
            sl5 <- lambda('U235')[2]
            E <- diag(c(st,sl5))^2
            R.err <- sqrt(J%*%E%*%t(J))
        } else {
            R.err <- abs(J*st)
        }
        out <- cbind(R,R.err)
        colnames(out) <- c('75','s[75]')
    }
    out
}
age_to_U235Pb207_ratio <- function(tt,st=0,d=diseq(),exterr=FALSE){
    ns <- length(tt)
    if (ns>1){
        out <- matrix(0,ns,2)
        colnames(out) <- c('57','s[57]')
        for (i in 1:ns){
            if (length(st)<ns) sti <- st[1]
            else sti <- st[i]
            out[i,] <- age_to_U235Pb207_ratio(tt[i],st=sti,d=d[i],exterr=exterr)
        }
    } else {
        l5 <- lambda('U235')[1]
        D <- mclean(tt=tt,d=d,exterr=exterr)
        R <- 1/D$Pb207U235
        if (exterr){
            J <- -c(D$dPb207U235dt,D$dPb207U235l35)/D$Pb207U235^2
            sl5 <- lambda('U235')[2]
            E <- diag(c(st,sl5))^2
            R.err <- sqrt(J%*%E%*%t(J))
        } else {
            J <- -D$dPb207U235dt/D$Pb207U235^2
            R.err <- abs(J*st)
        }
        out <- cbind(R,R.err)
        colnames(out) <- c('57','s[57]')
    }
    out
}
age_to_Pb206U238_ratio <- function(tt,st=0,d=diseq(),exterr=FALSE){
    ns <- length(tt)
    if (ns>1){
        out <- matrix(0,ns,2)
        colnames(out) <- c('68','s[68]')
        for (i in 1:ns){
            if (length(st)<ns) sti <- st[1]
            else sti <- st[i]
            out[i,] <- age_to_Pb206U238_ratio(tt[i],st=sti,d=d[i],exterr=exterr)
        }
    } else {
        l8 <- lambda('U238')[1]
        D <- mclean(tt=tt,d=d,exterr=exterr)
        R <- D$Pb206U238
        J <- D$dPb206U238dt
        if (exterr){
            J <- c(J,D$dPb206U238dl38)
            sl8 <- lambda('U238')[2]
            E <- diag(c(st,sl8))^2
            R.err <- sqrt(J%*%E%*%t(J))
        } else {
            R.err <- abs(J*st)
        }
        out <- cbind(R,R.err)
        colnames(out) <- c('68','s[68]')
    }
    out
}
age_to_U238Pb206_ratio <- function(tt,st=0,d=diseq(),exterr=FALSE){
    ns <- length(tt)
    if (ns>1){
        out <- matrix(0,ns,2)
        colnames(out) <- c('86','s[86]')
        for (i in 1:ns){
            if (length(st)<ns) sti <- st[1]
            else sti <- st[i]
            out[i,] <- age_to_U238Pb206_ratio(tt[i],st=sti,d=d[i],exterr=exterr)
        }
    } else {
        l8 <- lambda('U238')[1]
        tt <- check.zero.UPb(tt)
        D <- mclean(tt=tt,d=d,exterr=exterr)
        R <- 1/D$Pb206U238
        if (exterr){
            J <- -c(D$dPb206U238dt,D$dPb206U238l38)/D$Pb206U238^2
            sl8 <- lambda('U238')[2]
            E <- diag(c(st,sl8))^2
            R.err <- sqrt(J%*%E%*%t(J))
        } else {
            J <- -D$dPb206U238dt/D$Pb206U238^2
            R.err <- abs(J*st)
        }
        R.err <- abs(J*st)
        out <- cbind(R,R.err)
        colnames(out) <- c('86','s[86]')
    }
    out
}
age_to_Pb207Pb206_ratio <- function(tt,st=0,d=diseq(),exterr=FALSE){
    ns <- length(tt)
    if (ns>1){
        out <- matrix(0,ns,2)
        colnames(out) <- c('76','s[76]')
        for (i in 1:ns){
            if (length(st)<ns) sti <- st[1]
            else sti <- st[i]
            out[i,] <- age_to_Pb207Pb206_ratio(tt[i],st=sti,d=d[i],exterr=exterr)
        }
    } else {
        l5 <- lambda('U235')[1]
        l8 <- lambda('U238')[1]
        tt <- check.zero.UPb(tt)
        U <- iratio('U238U235')[1]
        D <- mclean(tt=tt,d=d,exterr=exterr)
        R <- (1/U)*D$Pb207U235/D$Pb206U238
        J <- (1/U)*(D$dPb207U235dt*D$Pb206U238 -
                    D$Pb207U235*D$dPb206U238dt)/D$Pb206U238^2
        if (exterr){
            sl5 <- lambda('U235')[2]
            sl8 <- lambda('U238')[2]
            d76dl35 <- (1/U)*D$dPb207U235dl35/D$Pb206U238
            d76dl38 <- -(1/U)*D$Pb207U235*D$dPb206U238dl38/D$Pb206U238^2
            J <- c(J,d76dl35,d76dl38)
            E <- diag(c(st,sl5,sl8))^2
            R.err <- sqrt(J%*%E%*%t(J))
        } else {
            R.err <- abs(J*st)
        }
        out <- cbind(R,R.err)
        colnames(out) <- c('76','s[76]')
    }
    out
}
age_to_Pb208Th232_ratio <- function(tt,st=0,exterr=FALSE){
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
        if (exterr){
            sl2 <- lambda('Th232')[2]
            J <- c(J,tt*exp(l2*tt))
            E <- diag(c(st,sl2))^2
            R.err <- sqrt(J%*%E%*%t(J))
        } else {
            R.err <- abs(J*st)
        }
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
    if (x$format %in% c(4,6)){
        out <- subset(x$x,select=labels)
    } else if (x$format == 5){
        Pb204U238 <- x$x[,'Pb204Pb206']/x$x[,'U238Pb206']
        errPb204U238 <-
            Pb204U238*sqrt( (x$x[,'errPb204Pb206']/x$x[,'Pb204Pb206'])^2 +
                            (x$x[,'errU238Pb206']/x$x[,'U238Pb206'])^2 )
        out <- cbind(Pb204U238,errPb204U238)
    } else {
        stop('No 204Pb measurements available!')
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
    labels <- c('Pb208Th232','errPb208Th232')
    if (x$format == 7){
        out <- x$x[,labels,drop=FALSE]
    } else if (x$format == 8){
        ns <- length(x)
        out <- matrix(0,ns,2)
        colnames(out) <- labels
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
    labels <- c('Pb208Pb206','errPb208Pb206')
    if (x$format == 7){
        ns <- length(x)
        out <- matrix(0,ns,2)
        colnames(out) <- labels
        out[,'Pb208Pb206'] <- x$x[,'Pb208Th232']*x$x[,'Th232U238']/x$x[,'Pb206U238']
        J1 <- -out[,'Pb208Pb206']/x$x[,'Pb206U238']
        J2 <- x$x[,'Th232U238']/x$x[,'Pb206U238']
        J3 <- x$x[,'Pb208Th232']/x$x[,'Pb206U238']
        E11 <- x$x[,'errPb206U238']^2
        E22 <- x$x[,'errPb208Th232']^2
        E33 <- x$x[,'Th232U238']^2
        E12 <- x$x[,'rhoXZ']*x$x[,'Pb206U238']*x$x[,'Pb208Th232']
        E13 <- x$x[,'rhoXW']*x$x[,'Pb206U238']*x$x[,'Th232U238']
        E23 <- x$x[,'rhoZW']*x$x[,'Pb208Th232']*x$x[,'Th232U238']
        out[,'errPb208Th232'] <- errorprop1x3(J1,J2,J3,E11,E22,E33,E12,E13,E23)
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
        t.init <- ifelse(x>-1,log(1+x)/l5,0)
        E <- matrix(0,3,3)
        J <- matrix(0,1,3)
        E[1,1] <- sx^2
        E[2:3,2:3] <- getEl('U235')
        if (d$equilibrium | d$PaU$option<1){
            t.75 <- t.init
            J[1,1] <- 1/(l5*(1+x))                       # dt/dx
            if (exterr & x>-1) J[1,2] <- log(1+x)/l5^2   # dt/dl5
        } else { # apply a disequilibrium correction
            t.75 <- stats::optimise(diseq.75.misfit,x=x,d=d,
                                    interval=t.init*c(.5,2))$minimum
            D <- mclean(tt=t.75,d=d,exterr=exterr)    
            J[1,1] <- 1/D$dPb207U235dt                 # dt/dx
            J[1,2] <- -D$dPb207U235dl35/D$dPb207U235dt # dt/dl35
            J[1,3] <- -D$dPb207U235dl31/D$dPb207U235dt # dt/dl31
        }
        st.75 <- sqrt(J %*% E %*% t(J))
        out <- c(t.75,st.75)
        names(out) <- c('t75','s[t75]')
    }
    out
}
get.Pb207U235.age.UPb <- function(x,i=1,exterr=FALSE,...){
    r75 <- get.Pb207U235.ratios(x)
    get.Pb207U235.age(r75[i,'Pb207U235'],r75[i,'errPb207U235'],
                      exterr=exterr,d=x$d[i])
}
get.Pb207U235.age.wetherill <- function(x,exterr=FALSE,...){
    i <- 'Pb207U235'
    r75 <- x$x[i]
    sr75 <- sqrt(x$cov[i,i])
    get.Pb207U235.age(r75,sr75,exterr=exterr,d=x$d,...)
}

get.Pb206U238.age <- function(x,...){ UseMethod("get.Pb206U238.age",x) }
get.Pb206U238.age.default <- function(x,sx=0,exterr=FALSE,d=diseq()){
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
        E <- matrix(0,5,5)
        E[1,1] <- sx^2
        E[2:5,2:5] <- getEl('U238')
        J <- matrix(0,1,5)
        l8 <- lambda('U238')[1]
        if (x>-1) t.init <- log(1+x)/l8 else t.init <- 0
        if (d$equilibrium){
            t.68 <- t.init
            J[1,1] <- 1/(l8*(1+x))                       # dt/dx
            if (exterr & x>-1) J[1,2] <- log(1+x)/l8^2   # dt/dl38
        } else { # apply a disequilibrium correction
            if (measured.disequilibrium(d)) tlim <- c(0,meas.diseq.maxt(d))
            else tlim <- c(0,4600)
            t.68 <- tryCatch({
                stats::optimise(diseq.68.misfit,interval=tlim,x=x,d=d)$minimum
            }, error = function(error_condition) {
                t.init
            })
            D <- mclean(tt=t.68,d=d,exterr=exterr)
            J[1,1] <- 1/D$dPb206U238dt
            J[1,2] <- -D$dPb206U238dl38/D$dPb206U238dt # dt/dl38
            J[1,3] <- -D$dPb206U238dl34/D$dPb206U238dt # dt/dl34
            J[1,4] <- -D$dPb206U238dl30/D$dPb206U238dt # dt/dl30
            J[1,5] <- -D$dPb206U238dl26/D$dPb206U238dt # dt/dl26
        }
        st.68 <- sqrt(J %*% E %*% t(J))
        out <- c(t.68,st.68)
        names(out) <- c('t68','s[t68]')
    }
    out
}
get.Pb206U238.age.UPb <- function(x,i=NULL,exterr=FALSE,...){
    r68 <- get.Pb206U238.ratios(x)
    if (is.null(i)){
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
get.Pb207Pb206.age.default <- function(x,sx=0,exterr=FALSE,d=diseq(),t.68=NULL,...){
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
        if (!d$equilibrium & !is.null(t.68)){
            midpoint <- stats::optimise(twslope,d=d,interval=interval)$minimum
            if (t.68<midpoint){
                interval[2] <- midpoint
            } else {
                interval[1] <- midpoint
            }
        }
        if (measured.disequilibrium(d)) interval[2] <- meas.diseq.maxt(d)
        t.76 <- stats::optimise(get.76.misfit,x=x,d=d,interval=interval)$minimum
        D <- mclean(tt=t.76,d=d,exterr=exterr)
        E <- matrix(0,9,9)
        E[1,1] <- sx^2
        E[2,2] <- iratio('U238U235')[2]^2
        E[3:9,3:9] <- getEl()
        J <- matrix(0,1,9)
        J[1,1] <- -1/D$dPb207Pb206dt                # dt/dx
        J[1,2] <- D$dPb207Pb206dU/D$dPb207Pb206dt   # dt/dU
        J[1,3] <- D$dPb207Pb206dl38/D$dPb207Pb206dt # dt/dl38
        J[1,4] <- D$dPb207Pb206dl35/D$dPb207Pb206dt # dt/dl35
        J[1,5] <- D$dPb207Pb206dl34/D$dPb207Pb206dt # dt/dl34
        J[1,7] <- D$dPb207Pb206dl31/D$dPb207Pb206dt # dt/dl31
        J[1,8] <- D$dPb207Pb206dl30/D$dPb207Pb206dt # dt/dl30
        J[1,9] <- D$dPb207Pb206dl26/D$dPb207Pb206dt # dt/dl26
        st.76 <- sqrt( J %*% E %*% t(J) )
        out <- c(t.76,st.76)
        names(out) <- c('t76','s[t76]')
    }
    out    
}
get.Pb207Pb206.age.UPb <- function(x,i=1,exterr=FALSE,...){
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
get.Pb208Th232.age.UPb <- function(x,i=NULL,exterr=FALSE,...){
    r82 <- get.Pb208Th232.ratios(x)
    if (is.null(i)){
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
UPb.age <- function(x,exterr=FALSE,i=NULL,conc=TRUE,omit4c=NULL,
                    discordance=discfilter(),common.Pb=0,...){
    if (discordance$option==0 | discordance$before) xd <- x
    else xd <- Pb0corr(x,option=common.Pb,omit4c=omit4c)
    if (common.Pb==0) X <- x
    else X <- Pb0corr(x,option=common.Pb,omit4c=omit4c)
    if (!is.null(i)){
        out <- UPb_age_helper(x=x,X=X,xd=xd,i=i,exterr=exterr,
                              conc=conc,discordance=discordance)
    } else {
        nn <- length(x)
        out <- NULL
        for (i in 1:nn){
            ti <- UPb_age_helper(x=x,X=X,xd=xd,i=i,exterr=exterr,
                                 conc=conc,discordance=discordance)
            out <- rbind(out,ti)
        }
    }
    out
}

# x = raw data
# X = common Pb corrected data (if common.Pb>0)
# xd = data to be used for concordia age calculation 
#      (raw if before==TRUE, common Pb corrected if before==FALSE)
UPb_age_helper <- function(x,X,xd,i=1,exterr=FALSE,
                           conc=TRUE,discordance=discfilter(),...){
    Xi <- subset(X,subset=((1:length(X))%in%i))
    labels <- c('t.75','s[t.75]','t.68','s[t.68]','t.76','s[t.76]')
    hasTh <- (x$format>6)
    if (hasTh) labels <- c(labels,'t.82','s[t.82]')
    if (conc) labels <- c(labels,'t.conc','s[t.conc]')
    tlabels <- labels
    if (discordance$option%in%c(1,'t',2,'r',3,'sk',4,'a',5,'c'))
        labels <- c(labels,'disc')
    if (discordance$option%in%c(6,'p'))
        labels <- c(labels,'p[conc]')
    t.75 <- get.Pb207U235.age(Xi,exterr=exterr)
    t.68 <- get.Pb206U238.age(Xi,exterr=exterr)
    t.76 <- get.Pb207Pb206.age(Xi,exterr=exterr,t.68=subset(t.68,select=1))
    out <- c(t.75,t.68,t.76)
    if (hasTh){
        t.82 <- get.Pb208Th232.age(Xi,exterr=exterr)
        out <- c(out,t.82)
    }
    if (conc){
        t.conc <- concordia.age(x=Xi,i=1,exterr=exterr)
        out <- c(out,t.conc$age)
    }
    if (discordance$option>0){
        xdi <- subset(xd,subset=((1:length(xd))%in%i))
    }
    if (discordance$option%in%c(1,'t',2,'r',3,'sk',4,'a',5,'c')){
        xi <- subset(x,subset=((1:length(x))%in%i))
        dif <- discordance(x=xi,X=xdi,option=discordance$option)
        out <- c(out,dif)
    }
    if (discordance$option%in%c(6,'p')){
        t.conc <- concordia.age(x=xdi,exterr=exterr)
        SS.concordance <-
            LL.concordia.age(pars=t.conc$age[1],cc=wetherill(xdi,i=1),
                             mswd=TRUE,exterr=exterr,d=xdi$d)
        p.value <- 1-stats::pchisq(SS.concordance,1)
        out <- c(out,p.value)
    }
    names(out) <- labels
    out
}

getEl <- function(parent){
    out <- matrix(0,7,7)
    out[1,1] <- lambda('U238')[2]^2
    out[2,2] <- lambda('U235')[2]^2
    out[3,3] <- (lambda('U234')[2]*1000)^2
    out[4,4] <- lambda('Th232')[2]^2
    out[5,5] <- (lambda('Pa231')[2]*1000)^2
    out[6,6] <- (lambda('Th230')[2]*1000)^2
    out[7,7] <- (lambda('Ra226')[2]*1000)^2
    rcnames <- c('U238','U235','U234','Th232','Pa231','Th230','Ra226')
    rownames(out) <- rcnames
    colnames(out) <- rcnames
    if (missing(parent)){
        # do nothing
    } else if (identical(parent,'U238')){
        rcnames <- c('U238','U234','Th230','Ra226')
    } else if (identical(parent,'U235')){
        rcnames <- c('U235','Pa231')
    } else if (identical(parent,'Th232')){
        rcnames <- c('Th232')
    } else {
        # do nothing
    }
    out[rcnames,rcnames]
}
