wetherill <- function(x,i=NA,exterr=FALSE){
    out <- list()
    if (x$format %in% c(1,3)){
        labels <- c('Pb207U235','Pb206U238')
        out$x <- x$x[i,labels]
        out$cov <- matrix(0,2,2)
        diag(out$cov) <-
            x$x[i,c('errPb207U235','errPb206U238')]^2
        out$cov <-
            cor2cov2(x$x[i,'errPb207U235'],x$x[i,'errPb206U238'],x$x[i,'rhoXY'])
    } else if (x$format == 2){
        labels <- c('Pb207U235','Pb206U238')
        U238U235 <- iratio('U238U235')[1]
        Pb207U235 <- U238U235*x$x[i,'Pb207Pb206']/x$x[i,'U238Pb206']
        Pb206U238 <- 1/x$x[i,'U238Pb206']
        out$x <- c(Pb207U235,Pb206U238)
        J <- matrix(0,2,3)
        E <- matrix(0,3,3)
        E[3,3] <- iratio('U238U235')[2]^2
        J[1,1] <- -Pb207U235/x$x[i,'U238Pb206']
        J[1,2] <- U238U235/x$x[i,'U238Pb206']
        if (exterr) J[1,3] <- x$x[i,'Pb207Pb206']/x$x[i,'U238Pb206']
        J[2,1] <- -1/x$x[i,'U238Pb206']^2
        E[1:2,1:2] <- cor2cov2(x$x[i,'errU238Pb206'],
                              x$x[i,'errPb207Pb206'],x$x[i,'rhoXY'])
        out$cov <- J %*% E %*% t(J)
    } else if (x$format == 4){
        labels <- c('Pb207U235','Pb206U238','Pb204U238')
        out$x <- x$x[i,labels]
        out$cov <- matrix(0,3,3)
        diag(out$cov) <-
            x$x[i,c('errPb207U235','errPb206U238','errPb204U238')]^2
        out$cov <-
            cor2cov3(x$x[i,'errPb207U235'],x$x[i,'errPb206U238'],
                     x$x[i,'errPb204U238'],x$x[i,'rhoXY'],
                     x$x[i,'rhoXZ'],x$x[i,'rhoYZ'])
    } else if (x$format == 5){
        labels <- c('Pb207U235','Pb206U238','Pb204U238')
        U238U235 <- iratio('U238U235')[1]
        Pb207U235 <- U238U235*x$x[i,'Pb207Pb206']/x$x[i,'U238Pb206']
        Pb206U238 <- 1/x$x[i,'U238Pb206']
        Pb204U238 <- x$x[i,'Pb204Pb206']/x$x[i,'U238Pb206']
        out$x <- c(Pb207U235,Pb206U238,Pb204U238)
        J <- matrix(0,3,4)
        E <- matrix(0,4,4)
        E[4,4] <- iratio('U238U235')[2]^2
        J[1,1] <- -Pb207U235/x$x[i,'U238Pb206']
        J[1,2] <- U238U235/x$x[i,'U238Pb206']
        if (exterr) J[1,4] <- x$x[i,'Pb207Pb206']/x$x[i,'U238Pb206']
        J[2,1] <- -1/x$x[i,'U238Pb206']^2
        J[3,1] <- Pb204U238/x$x[i,'U238Pb206']
        J[3,3] <- 1/x$x[i,'U238Pb206']
        E[1:3,1:3] <- cor2cov3(x$x[i,'errU238Pb206'],x$x[i,'errPb207Pb206'],
                               x$x[i,'errPb204Pb206'],x$x[i,'rhoXY'],
                               x$x[i,'rhoXZ'],x$x[i,'rhoYZ'])
        out$cov <- J %*% E %*% t(J)        
    } else if (x$format == 6){
        labels <- c('Pb207U235','Pb206U238','Pb204U238')
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
    }
    names(out$x) <- labels
    colnames(out$cov) <- labels
    rownames(out$cov) <- labels
    class(out) <- "wetherill"
    out
}
tera.wasserburg <- function(x,i,exterr=FALSE){
    out <- list()
    if (x$format==1){
        labels <- c('U238Pb206','Pb207Pb206')
        U238U235 <- iratio('U238U235')[1]
        U238Pb206 <- 1/x$x[i,'Pb206U238']
        Pb207Pb206 <- x$x[i,'Pb207U235']/(x$x[i,'Pb206U238']*U238U235)
        J <- matrix(0,2,3)
        E <- matrix(0,3,3)
        E[3,3] <- iratio('U238U235')[2]^2
        J[1,2] <- -1/x$x[i,'Pb206U238']^2
        J[2,1] <- 1/(x$x[i,'Pb206U238']*U238U235)
        J[2,2] <- -Pb207Pb206/x$x[i,'Pb206U238']
        if (exterr) J[2,3] <- -Pb207Pb206/U238U235
        E[1:2,1:2] <- cor2cov2(x$x[i,'errPb207U235'],
                              x$x[i,'errPb206U238'],x$x[i,'rhoXY'])
        out$x <- c(U238Pb206,Pb207Pb206)
        out$cov <- J %*% E %*% t(J)
    } else if (x$format == 2){
        labels <- c('U238Pb206','Pb207Pb206')
        out$x <- x$x[i,c('U238Pb206','Pb207Pb206')]
        out$cov <- matrix(0,2,2)
        diag(out$cov) <- x$x[i,c('errU238Pb206','errPb207Pb206')]^2
        out$cov <-
            cor2cov2(x$x[i,'errU238Pb206'],x$x[i,'errPb207Pb206'],x$x[i,'rhoXY'])
    } else if (x$format==3){
        labels <- c('U238Pb206','Pb207Pb206')
        U238Pb206 <- 1/x$x[i,'Pb206U238']
        out$x <- c(U238Pb206,x$x[i,'Pb207Pb206'])
        J <- matrix(0,2,2)
        E <- matrix(0,2,2)
        E <- cor2cov2(x$x[i,'errPb206U238'],x$x[i,'errPb207Pb206'],x$x[i,'rhoYZ'])
        J[1,1] <- -U238Pb206^2
        J[2,2] <- 1
        out$cov <- J %*% E %*% t(J)
    } else if (x$format == 4){
        labels <- c('U238Pb206','Pb207Pb206','Pb204Pb206')
        U238U235 <- iratio('U238U235')[1]
        U238Pb206 <- 1/x$x[i,'Pb206U238']
        Pb207Pb206 <- x$x[i,'Pb207U235']/(x$x[i,'Pb206U238']*U238U235)
        Pb204Pb206 <- x$x[i,'Pb204U238']/x$x[i,'Pb206U238']
        J <- matrix(0,3,4)
        E <- matrix(0,4,4)
        E[4,4] <- iratio('U238U235')[2]^2
        J[1,2] <- -1/x$x[i,'Pb206U238']^2
        J[2,1] <- 1/(x$x[i,'Pb206U238']*U238U235)
        J[2,2] <- -Pb207Pb206/x$x[i,'Pb206U238']
        if (exterr) J[2,4] <- -Pb207Pb206/U238U235
        J[3,2] <- -Pb204Pb206/x$x[i,'Pb206U238']
        J[3,3] <- 1/x$x[i,'Pb206U238']
        E[1:3,1:3] <- cor2cov3(x$x[i,'errPb207U235'],x$x[i,'errPb206U238'],
                               x$x[i,'errPb204U238'],x$x[i,'rhoXY'],
                               x$x[i,'rhoXZ'],x$x[i,'rhoXZ'])
        out$x <- c(U238Pb206,Pb207Pb206,Pb204Pb206)
        out$cov <- J %*% E %*% t(J)
    } else if (x$format == 5){
        labels <- c('U238Pb206','Pb207Pb206','Pb204Pb206')
        out$x <- x$x[i,c('U238Pb206','Pb207Pb206','Pb204Pb206')]
        out$cov <- matrix(0,3,3)
        diag(out$cov) <- x$x[i,c('errU238Pb206','errPb207Pb206','errPb204Pb206')]^2
        out$cov <-
            cor2cov3(x$x[i,'errU238Pb206'],x$x[i,'errPb207Pb206'],x$x[i,'errPb204Pb206'],
                     x$x[i,'rhoXY'],x$x[i,'rhoXZ'],x$x[i,'rhoYZ'])
    } else if (x$format == 6){
        labels <- c('U238Pb206','Pb207Pb206','Pb204Pb206')
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
    }
    names(out$x) <- labels
    colnames(out$cov) <- labels
    rownames(out$cov) <- labels
    class(out) <- "terawasserburg"
    out
}

# convert data to a 5-column table for concordia analysis
flat.UPb.table <- function(x,wetherill=TRUE){
    ns <- length(x)
    out <- matrix(0,ns,5)
    for (i in 1:ns){
        if (wetherill) xi <- wetherill(x,i)
        else xi <- tera.wasserburg(x,i)
        out[i,1] <- xi$x[1]
        out[i,3] <- xi$x[2]
        out[i,2] <- sqrt(xi$cov[1,1])
        out[i,4] <- sqrt(xi$cov[2,2])
        out[i,5] <- xi$cov[1,2]/(out[i,2]*out[i,4])
    }
    if (wetherill){
        colnames(out) <- c('Pb207U235','errPb207U235',
                           'Pb206U238','errPb206U238','rhoXY')
        class(out) <- 'wetherill'
    } else {
        colnames(out) <- c('U238Pb206','errU238Pb206',
                           'Pb207Pb206','errPb207Pb206','rhoXY')
        class(out) <- 'terawasserburg'
    }
    out
}

age_to_wetherill_ratios <- function(tt,st=0,exterr=FALSE){
    out <- list()
    labels <- c('Pb207U235','Pb206U238')
    l8 <- settings('lambda','U238')[1]
    l5 <- settings('lambda','U235')[1]
    out$x <- c(exp(l5*tt)-1,exp(l8*tt)-1)
    E <- matrix(0,3,3)
    diag(E) <- c(st,lambda('U235')[2],lambda('U238')[2])^2
    J <- matrix(0,2,3)
    J[1,1] <- l5*exp(l5*tt)
    J[2,1] <- l8*exp(l8*tt)
    if (exterr){
        J[1,2] <- tt*exp(l5*tt)
        J[2,3] <- tt*exp(l8*tt)
    }
    out$cov <- J %*% E %*% t(J)
    names(out$x) <- labels
    rownames(out$cov) <- labels
    colnames(out$cov) <- labels
    out
}
age_to_terawasserburg_ratios <- function(tt,st=0,exterr=FALSE){
    out <- list()
    labels <- c('U238Pb206','Pb207Pb206')
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    tt <- check.zero.UPb(tt)
    R <- iratio('U238U235')[1]
    Pb207U235 <- exp(l5*tt)-1
    Pb206U238 <- exp(l8*tt)-1
    U238Pb206 <- 1/Pb206U238
    Pb207Pb206 <- (1/R)*(exp(l5*tt)-1)/(exp(l8*tt)-1)
    out$x <- c(U238Pb206,Pb207Pb206)
    E <- matrix(0,4,4)
    diag(E) <- c(st,lambda('U235')[2],lambda('U238')[2],iratio('U238U235')[2])^2
    J <- matrix(0,2,4)
    J[1,1] <- -l8*U238Pb206^2
    J[2,1] <- (1/R)*(l5*exp(l5*tt)*Pb206U238-l8*exp(l8*tt)*Pb207U235)*U238Pb206^2
    if (exterr){
        J[1,3] <- -tt*exp(l8*tt)*U238Pb206^2
        J[2,2] <- tt*exp(l5*tt)*U238Pb206/R
        J[2,3] <- tt*Pb207U235*U238Pb206/R
        J[2,4] <- -Pb207Pb206/R
    }
    out$cov <- J %*% E %*% t(J)
    names(out$x) <- labels
    rownames(out$cov) <- labels
    colnames(out$cov) <- labels
    out
}

age_to_Pb207U235_ratio <- function(tt,st=0){
    l5 <- lambda('U235')[1]
    R <- exp(l5*tt)-1
    J <- l5*exp(l5*tt)
    R.err <- J*st
    out <- cbind(R,R.err)
    colnames(out) <- c('75','s[75]')
    out
}
age_to_Pb206U238_ratio <- function(tt,st=0){
    l8 <- lambda('U238')[1]
    R <- exp(l8*tt)-1
    J <- l8*exp(l8*tt)
    R.err <- J*st
    out <- cbind(R,R.err)
    colnames(out) <- c('68','s[68]')
    out
}
age_to_U238Pb206_ratio <- function(tt,st=0){
    l8 <- lambda('U238')[1]
    tt <- check.zero.UPb(tt)
    R <- 1/(exp(l8*tt)-1)
    J <- -l8*exp(l8*tt)/(exp(l8*tt)-1)^2
    R.err <- abs(J*st)
    out <- cbind(R,R.err)
    colnames(out) <- c('86','s[86]')
    out
}
age_to_Pb207Pb206_ratio <- function(tt,st=0){
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    tt <- check.zero.UPb(tt)
    U <- iratio('U238U235')[1]
    R <- (1/U)*(exp(l5*tt)-1)/(exp(l8*tt)-1)
    N <- l5*exp(l5*tt)*(exp(l8*tt)-1) - l8*exp(l8*tt)*(exp(l5*tt)-1)
    D <- (exp(l8*tt)-1)^2
    J <- (1/U)*(N/D)
    R.err <- abs(J*st)
    out <- cbind(R,R.err)
    colnames(out) <- c('76','s[76]')
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
    if (x$format < 4) stop('No 204Pb measurements available!')
    else if (x$format %in% c(4,6)){
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
    if (x$format %in% c(1,3,4,6)){
        out <- subset(x$x,select=labels)
    } else if (x$format %in% c(2,5)){
        R <- iratio('U238U235')[1]
        sR <- iratio('U238U235')[2]
        X <- x$x[,'U238Pb206']
        sX <- x$x[,'errU238Pb206']
        Y <- x$x[,'Pb207Pb206']
        sY <- x$x[,'errPb207Pb206']
        rhoXY <- x$x[,'rhoXY']
        covXY <- rhoXY*sX*sY
        out[,'Pb207U235'] <- R*Y/X
        relerr2 <- (sX/X)^2 -2*covXY/(X*Y) + (sY/Y)^2
        if (exterr) relerr2 <- relerr2 + (sR/R)^2
        out[,'errPb207U235'] <- sqrt(relerr2)*out[,'Pb207U235']
    }
    out
}
get.Pb206U238.ratios <- function(x){
    ns <- length(x)
    out <- matrix(0,ns,2)
    labels <- c('Pb206U238','errPb206U238')
    colnames(out) <- labels
    if (x$format %in% c(1,3,4,6)){
        out <- subset(x$x,select=labels)
    } else if (x$format %in% c(2,5)){
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
    if (x$format %in% c(1,3,4,6)){
        out[,'U238Pb206'] <- 1/x$x[,'Pb206U238']
        out[,'errU238Pb206'] <- out[,'U238Pb206']*
            x$x[,'errPb206U238']/x$x[,'Pb206U238']
    } else if (x$format %in% c(2,5)){
        out <- subset(x$x,select=labels)
    }
    out
}
get.Pb207Pb206.ratios <- function(x,exterr=FALSE){
    ns <- length(x)
    out <- matrix(0,ns,2)
    labels <- c('Pb207Pb206','errPb207Pb206')
    colnames(out) <- labels
    if (x$format %in% c(1,4)){
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
    } else if (x$format %in% c(2,3,5,6)){
        out <- subset(x$x,select=labels)
    }
    out
}

get.Pb207U235.age <- function(x,...){ UseMethod("get.Pb207U235.age",x) }
get.Pb207U235.age.default <- function(x,sx=0,exterr=TRUE,...){
    l5 <- lambda('U235')[1]
    sl5 <- lambda('U235')[2]
    t.75 <- log(1+x)/l5
    J <- matrix(0,1,2)
    J[1,1] <- 1/(l5*(1+x))
    if (exterr) J[1,2] <- log(1+x)/l5^2
    E <- matrix(0,2,2)
    E[1,1] <- sx^2
    E[2,2] <- sl5^2
    st.75 <- sqrt(J %*% E %*% t(J))
    out <- c(t.75,st.75)
    names(out) <- c('t75','s[t75]')
    out
}
get.Pb207U235.age.UPb <- function(x,i,exterr=TRUE,...){
    r75 <- get.Pb207U235.ratios(x)
    get.Pb207U235.age(r75[i,'Pb207U235'],r75[i,'errPb207U235'],exterr=exterr)
}
get.Pb207U235.age.wetherill <- function(x,exterr=TRUE,...){
    i <- 'Pb207U235'
    r75 <- x$x[i]
    sr75 <- sqrt(x$cov[i,i])
    get.Pb207U235.age(r75,sr75,exterr=exterr,...)
}

get.Pb206U238.age <- function(x,...){ UseMethod("get.Pb206U238.age",x) }
get.Pb206U238.age.default <- function(x,sx=0,exterr=TRUE,...){
    l8 <- lambda('U238')[1]
    sl8 <- lambda('U238')[2]
    t.68 <- log(1+x)/l8
    if (length(x)>1){
        out <- t.68
    } else {
        J <- matrix(0,1,2)
        J[1,1] <- 1/(l8*(1+x))
        if (exterr) J[1,2] <- log(1+x)/l8^2
        E <- matrix(0,2,2)
        E[1,1] <- sx^2
        E[2,2] <- sl8^2
        st.68 <- sqrt(J %*% E %*% t(J))
        out <- c(t.68,st.68)
        names(out) <- c('t68','s[t68]')
    }
    out
}
get.Pb206U238.age.UPb <- function(x,i=NA,exterr=TRUE,...){
    r68 <- get.Pb206U238.ratios(x)
    if (is.na(i)){
        ns <- length(x)
        out <- matrix(0,ns,2)
        for (j in 1:ns){
            out[j,] <- get.Pb206U238.age.UPb(x,i=j,exterr=exterr,...)
        }
    } else {
        out <- get.Pb206U238.age(r68[i,'Pb206U238'],r68[i,'errPb206U238'],exterr=exterr)
    }
    out
}
get.Pb206U238.age.wetherill <- function(x,exterr=TRUE,...){
    i <- 'Pb206U238'
    r68 <- x$x[i]
    sr68 <- sqrt(x$cov[i,i])
    get.Pb206U238.age(r68,sr68,exterr=exterr,...)
}
get.Pb206U238.age.terawasserburg <- function(x,exterr=TRUE,...){
    i <- 'U238Pb206'
    r86 <- x$x[i]
    r68 <- 1/r86
    sr68 <- sqrt(x$cov[i,i])/r86
    get.Pb206U238.age(r68,sr68,exterr=exterr,...)
}

get.Pb207Pb206.age <- function(x,...){ UseMethod("get.Pb207Pb206.age",x) }
get.Pb207Pb206.age.default <- function(x,sx=0,exterr=TRUE,...){
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    sl5 <- lambda('U235')[2]
    sl8 <- lambda('U238')[2]
    R <- iratio('U238U235')[1]
    sR <- iratio('U238U235')[2]
    Pb207Pb206.misfit <- function(tt,x) {
        (age_to_Pb207Pb206_ratio(tt)[,'76']-x)^2
    }
    if (is.na(x)){
        t.76 <- NA
    } else {
        fit <- stats::optimize(Pb207Pb206.misfit,c(0,1e4),x=x)
        t.76 <- fit$minimum
    }
    J <- matrix(0,1,4)
    J[1,1] <- -(-1)/dD76dt(t.76,l5,l8,R)                      # dt/d76
    if (exterr) {
        J[1,2] <- -dD76dl5(t.76,l5,l8,R)/dD76dt(t.76,l5,l8,R) # dt/dl5
        J[1,3] <- -dD76dl8(t.76,l5,l8,R)/dD76dt(t.76,l5,l8,R) # dt/dl8
        J[1,4] <- -dD76dR(t.76,l5,l8,R)/dD76dt(t.76,l5,l8,R)  # dt/dR
    }
    E <- matrix(0,4,4)
    E[1,1] <- sx^2
    E[2,2] <- sl5^2
    E[3,3] <- sl8^2
    E[4,4] <- sR^2
    st.76 <- sqrt( J %*% E %*% t(J) )
    out <- c(t.76,st.76)
    names(out) <- c('t76','s[t76]')
    out    
}
get.Pb207Pb206.age.UPb <- function(x,i,exterr=TRUE,...){
    r76 <- get.Pb207Pb206.ratios(x)
    get.Pb207Pb206.age(r76[i,'Pb207Pb206'],
                       r76[i,'errPb207Pb206'],exterr=exterr)
}
get.Pb207Pb206.age.wetherill <- function(x,exterr=TRUE,...){
    U238U235 <- iratio('U238U235')[1]
    r76 <- x$x['Pb207U235']/(x$x['Pb206U238']*U238U235)
    J <- matrix(0,1,2)
    E <- x$cov
    J[1,1] <- 1/(x$x['Pb206U238']*U238U235)
    J[1,2] <- -x$x['Pb207U235']/(U238U235*x$x['Pb206U238']^2)
    sr76 <- J %*% E %*% t(J)
    get.Pb207Pb206.age(r76,sr76,exterr=exterr,...)
}
get.Pb207Pb206.age.terawasserburg <- function(x,exterr=TRUE,...){
    r76 <- x$x['Pb207Pb206']
    sr76 <- x$cov['Pb207Pb206','Pb207Pb206']
    get.Pb207Pb206.age(r76,sr76,exterr=exterr,...)
}

# x is an object of class \code{UPb}
# returns a matrix of 7/5, 6/8, 7/6
# and concordia ages and their uncertainties.
UPb.age <- function(x,exterr=TRUE,i=NA,sigdig=NA,conc=TRUE,show.p=FALSE){
    labels <- c('t.75','s[t.75]','t.68','s[t.68]','t.76','s[t.76]')
    if (conc) labels <- c(labels,'t.conc','s[t.conc]')
    if (conc & show.p) labels <- c(labels,'p[conc]')
    if (!is.na(i)){
        t.75 <- get.Pb207U235.age(x,i,exterr=exterr)
        t.68 <- get.Pb206U238.age(x,i,exterr=exterr)
        t.76 <- get.Pb207Pb206.age(x,i,exterr=exterr)
        t.75.out <- roundit(t.75[1],t.75[2],sigdig=sigdig)
        t.68.out <- roundit(t.68[1],t.68[2],sigdig=sigdig)
        t.76.out <- roundit(t.76[1],t.76[2],sigdig=sigdig)
        out <- c(t.75.out,t.68.out,t.76.out)
        if (conc){
            t.conc <- concordia.age(x,i,exterr=exterr)
            t.conc.out <- roundit(t.conc$age[1],t.conc$age[2],sigdig=sigdig)
            SS.concordance <-
                LL.concordia.age(tt=t.conc$age[1],ccw=wetherill(x,i),
                                 mswd=TRUE,exterr=exterr)
            out <- c(out,t.conc.out)
        }
        if (conc & show.p){
            p.value <- 1-stats::pchisq(SS.concordance,1)
            if (!is.na(sigdig)) p.value <- signif(p.value,sigdig)
            out <- c(out,p.value)
        }
        names(out) <- labels
    } else {
        nn <- nrow(x$x)
        out <- matrix(0,nn,length(labels))
        for (i in 1:nn){
            out[i,] <- UPb.age(x,i=i,exterr=exterr,sigdig=sigdig,
                               conc=conc,show.p=show.p)
        }
        colnames(out) <- labels
    }
    out
}

dD76dt <- function(t.76,l5,l8,R){
    el8t1 <- exp(l8*t.76)-1
    el5t1 <- exp(l5*t.76)-1
    num <- el8t1*l5*exp(l5*t.76)-el5t1*l8*exp(l8*t.76)
    den <- R*el8t1^2
    num/den
}

dD76dl5 <- function(t.76,l5,l8,R){
    el8t1 <- exp(l8*t.76)-1
    t.76*exp(l5*t.76)/(R*el8t1)
}

dD76dl8 <- function(t.76,l5,l8,R){
    el8t1 <- exp(l8*t.76)-1
    el5t1 <- exp(l5*t.76)-1
    -t.76*exp(l8*t.76)*el5t1/(R*el8t1^2)
}

dD76dR <- function(t.76,l5,l8,R){
    el5t1 <- exp(l5*t.76)-1
    el8t1 <- exp(l8*t.76)-1
    -el5t1/(el8t1*R^2)
}

filter.UPb.ages <- function(x,type=4,cutoff.76=1100,
                            cutoff.disc=c(-15,5),exterr=TRUE){
    tt <- UPb.age(x,exterr=exterr,conc=(type==5))
    do.76 <- tt[,'t.68'] > cutoff.76
    if (any(is.na(cutoff.disc))) {
        is.concordant <- rep(TRUE,nrow(x))
    } else {
        disc.75.68 <- 100*(1-tt[,'t.75']/tt[,'t.68'])
        disc.68.76 <- 100*(1-tt[,'t.68']/tt[,'t.76'])
        is.concordant <- (disc.75.68>cutoff.disc[1] & disc.75.68<cutoff.disc[2]) |
                         (disc.68.76>cutoff.disc[1] & disc.68.76<cutoff.disc[2])
        if (!any(is.concordant))
            stop(paste0('There are no concordant grains in this sample.',
                        'Try adjusting the discordance limits OR ',
                        'apply a common-Pb correction.'))
    }
    out <- matrix(NA,length(x),2)
    if (type==1){
        out[is.concordant,] <- tt[is.concordant,c('t.75','s[t.75]')]
    } else if (type==2){
        out[is.concordant,] <- tt[is.concordant,c('t.68','s[t.68]')]
    } else if (type==3){
        out[is.concordant,] <- tt[is.concordant,c('t.76','s[t.76]')]
    } else if (type==4){
        i.76 <- as.vector(which(do.76 & is.concordant))
        i.68 <- as.vector(which(!do.76 & is.concordant))
        out[i.76,] <- tt[i.76,c('t.76','s[t.76]')]
        out[i.68,] <- tt[i.68,c('t.68','s[t.68]')]
    } else if (type==5){
        out[is.concordant,] <- tt[is.concordant,c('t.conc','s[t.conc]')]
    }
    colnames(out) <- c('t','s[t]')
    out
}
