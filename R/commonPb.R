# option = 1: Stacey-Kramers
# option = 2: isochron
# option = 3: assumed Pb-composition
common.Pb.correction <- function(x,option=1){
    ns <- length(x)
    if (option == 1)
        out <- common.Pb.stacey.kramers(x)
    else if (option == 2)
        out <- common.Pb.isochron(x)
    else if (option == 3)
        out <- common.Pb.nominal(x)
    else out <- x
    out$x.raw <- x$x
    out
}

common.Pb.isochron <- function(x){
    fit <- ludwig(x)
    if (x$format<4){
        rr <- age_to_terawasserburg_ratios(fit$par[1],d=x$d)$x
        slope <- (rr['Pb207Pb206']-fit$par['76i'])/rr['U238Pb206']
        y0 <- get.Pb207Pb206.ratios(x)[,1] -
            slope*get.U238Pb206.ratios(x)[,1]
        out <- Pb.correction.without.204(x,i76=y0)
    } else {
        w <- wendt(fit$par['t'],d=x$d)
        r68i <- age_to_Pb206U238_ratio(fit$par['t'],d=x$d)[1]
        y86 <- data2york(x,option=3) # 4/6 vs. 8/6
        r86 <- y86[,'X'] + y86[,'Y']*fit$par['64i']/r68i
        r75i <- age_to_Pb207U235_ratio(fit$par['t'],d=x$d)[1]
        y57 <- data2york(x,option=4) # 4/7 vs. 5/7
        r57 <- y57[,'X'] + y57[,'Y']*fit$par['74i']/r75i
        out <- x
        ns <- length(x)
        U <- iratio('U238U235')[1]
        J1 <- matrix(0,2,3)
        J1[1,1] <- 1                   # dr86d86
        J1[1,3] <- fit$par['64i']/r68i # dr86d46
        J1[2,3] <- fit$par['74i']/r75i # dr57d46
        J2 <- matrix(0,2,2)
        out$x.raw <- out$x
        out$x <- matrix(0,ns,5)
        colnames(out$x) <- c('U238Pb206','errU238Pb206',
                             'Pb207Pb206','errPb207Pb206','rhoXY')
        out$format <- 2
        for (i in 1:ns){
            tw <- tera.wasserburg(x,i=i)
            E1 <- tw$cov
            J1[2,1] <- 1/(tw$x['Pb207Pb206']*U)   # dr57d86
            J1[2,2] <- -r57[i]/tw$x['Pb207Pb206'] # dr57d57
            E2 <- J1 %*% E1 %*% t(J1)
            out$x[i,'U238Pb206'] <- r86[i]
            out$x[i,'Pb207Pb206'] <- r86[i]/(r57[i]*U)
            J2[1,1] <- 1                          # dr86dr86
            J2[2,1] <- 1/(r57[i]*U)               # dr76dr86
            J2[2,2] <- -r86[i]/(U*r57[i]^2)       # dr76dr57
            E3 <- J2 %*% E2 %*% t(J2)
            out$x[i,'errU238Pb206'] <- sqrt(E3[1,1])
            out$x[i,'errPb207Pb206'] <- sqrt(E3[2,2])
            out$x[i,'rhoXY'] <- stats::cov2cor(E3)[1,2]
        }
    }
    out
}

common.Pb.stacey.kramers <- function(x){
    if (x$format < 4){
        tinit <- 1000 # initial guess
        for (i in 1:5){
            i6474 <- stacey.kramers(tinit)
            out <- Pb.correction.without.204(x,i6474[2]/i6474[1])
            tinit <- get.Pb206U238.age(out)[,1]
        }
    } else {
        ns <- length(x)
        out <- x
        out$x.raw <- x$x
        out$x <- matrix(0,ns,5)
        for (i in 1:ns){
            fit <- stats::optimise(SS.SK,interval=c(0,5000),x=x,i=i)
            ccw <- sk2w(x,i,tt=fit$minimum)
            out$x[i,c(1,3)] <- ccw$x
            out$x[i,c(2,4)] <- sqrt(diag(ccw$cov))
            out$x[i,5] <- stats::cov2cor(ccw$cov)[1,2]
        }
        colnames(out$x) <- c('Pb207U235','errPb207U235',
                             'Pb206U238','errPb206U238','rhoXY')
        out$format <- 1
    }
    out
}
SS.SK <- function(tt,x,i){
    ccw <- sk2w(x,i,tt)
    LL.concordia.age(tt,ccw,mswd=FALSE,exterr=FALSE,d=x$d)
}
sk2w <- function(x,i,tt){
    wi <- wetherill(x,i=i)
    i6474 <- stacey.kramers(tt)
    i64 <- i6474[1]
    i74 <- i6474[2]
    U <- iratio('U238U235')[1]
    out <- list(x=rep(0,2),cov=matrix(0,2,2))
    out$x[1] <- wi$x['Pb207U235'] - i74*wi$x['Pb204U238']*U
    out$x[2] <- wi$x['Pb206U238'] - i64*wi$x['Pb204U238']
    J <- matrix(0,2,3)
    J[1,1] <- 1
    J[1,3] <- -i74*U
    J[2,2] <- 1
    J[2,3] <- -i64
    out$cov <- J %*% wi$cov %*% t(J)
    out
}

common.Pb.nominal <- function(x){
    if (x$format < 4){
        i76 <- settings('iratio','Pb207Pb206')[1]
        out <- Pb.correction.without.204(x,i76)
    } else {
        i64 <- settings('iratio','Pb206Pb204')[1]
        i74 <- settings('iratio','Pb207Pb204')[1]
        out <- Pb.correction.with.204(x,i64,i74)
    }
    out
}

Pb.correction.without.204 <- function(x,i76){
    ns <- length(x)
    ni <- length(i76)
    out <- x
    m76 <- get.Pb207Pb206.ratios(x)[,1]
    m86 <- get.U238Pb206.ratios(x)[,1]
    tint <- rep(0,ns) # intercept age
    i76i <- i76[1]
    for (i in 1:ns){
        if (ni>1) i76i <- i76[i]
        tint[i] <- project.concordia(m76[i],m86[i],i76i,d=x$d)
    }
    if (x$format == 1){
        out$x[,'Pb207U235'] <- age_to_Pb207U235_ratio(tint,d=x$d)[,'75']
        out$x[,'errPb207U235'] <- x$x[,'errPb207U235']
        out$x[,'Pb206U238'] <- age_to_Pb206U238_ratio(tint,d=x$d)[,'68']
        out$x[,'errPb206U238'] <- x$x[,'errPb206U238']
        out$x[,'rhoXY'] <- x$x[,'rhoXY']        
    } else if (x$format == 2){
        out$x[,'U238Pb206'] <- age_to_U238Pb206_ratio(tint,d=x$d)[,'86']
        out$x[,'errU238Pb206'] <- x$x[,'errU238Pb206']
        out$x[,'Pb207Pb206'] <- age_to_Pb207Pb206_ratio(tint,d=x$d)[,'76']
        out$x[,'errPb207Pb206'] <- x$x[,'errPb207Pb206']
        out$x[,'rhoXY'] <- x$x[,'rhoXY']
    } else if (x$format == 3){
        out$x[,'Pb207Pb206'] <- age_to_Pb207Pb206_ratio(tint,d=x$d)[,'76']
        out$x[,'errPb207Pb206'] <- x$x[,'errPb207Pb206']
        out$x[,'Pb207U235'] <- age_to_Pb207U235_ratio(tint,d=x$d)[,'75']
        out$x[,'errPb207U235'] <- x$x[,'errPb207U235']
        out$x[,'Pb206U238'] <- age_to_Pb206U238_ratio(tint,d=x$d)[,'68']
        out$x[,'errPb206U238'] <- x$x[,'errPb206U238']
    }
    out
}
Pb.correction.with.204 <- function(x,i64,i74){
    ns <- length(x)
    out <- x
    U238U235 <- settings('iratio','U238U235')[1]
    if (x$format == 4){
        out$x <- matrix(0,ns,5)
        colnames(out$x) <-
            c('Pb207U235','errPb207U235','Pb206U238','errPb206U238','rhoXY')
        out$x[,'Pb207U235'] <- x$x[,'Pb207U235'] - i74*x$x[,'Pb204U238']*U238U235
        out$x[,'errPb207U235'] <- x$x[,'errPb207U235']
        out$x[,'Pb206U238'] <- x$x[,'Pb206U238'] - i64*x$x[,'Pb204U238']
        out$x[,'errPb206U238'] <- x$x[,'errPb206U238']
        out$x[,'rhoXY'] <- x$x[,'rhoXY']
        out$format <- 1
    } else if (x$format == 5){
        out$x <- matrix(0,ns,5)
        colnames(out$x) <-
            c('U238Pb206','errU238Pb206','Pb207Pb206','errPb207Pb206','rhoXY')
        out$x[,'U238Pb206'] <- x$x[,'U238Pb206']/(1 - i64*x$x[,'Pb204Pb206'])
        out$x[,'errU238Pb206'] <- x$x[,'errU238Pb206']
        out$x[,'Pb207Pb206'] <-
            (x$x[,'Pb207Pb206'] - x$x[,'Pb204Pb206']*i74)/(1 - x$x[,'Pb204Pb206']*i64)
        out$x[,'errPb207Pb206'] <- x$x[,'errPb207Pb206']
        out$x[,'rhoXY'] <- x$x[,'rhoXY']
        out$format <- 2
    } else if (x$format == 6){
        out$x <- matrix(0,ns,9)
        colnames(out$x) <- c('Pb207Pb206','errPb207Pb206','Pb207U235','errPb207U235',
                             'Pb206U238','errPb206U238','rhoXY','rhoXZ','rhoYZ')
        out$x[,'Pb207Pb206'] <-
            (x$x[,'Pb207Pb206'] - x$x[,'Pb204Pb206']*i74)/(1 - x$x[,'Pb204Pb206']*i64)
        out$x[,'errPb207Pb206'] <- x$x[,'errPb207Pb206']
        out$x[,'Pb207U235'] <- x$x[,'Pb207U235'] - i74*x$x[,'Pb204U238']*U238U235
        out$x[,'errPb207U235'] <- x$x[,'errPb207U235']
        out$x[,'Pb206U238'] <- x$x[,'Pb206U238'] - i64*x$x[,'Pb204U238']
        out$x[,'errPb206U238'] <- x$x[,'errPb206U238']
        out$format <- 3
    }
    out
}

stacey.kramers <- function(tt,inverse=FALSE){
    nt <- length(tt)
    sk.206.204 <- rep(0,nt)
    sk.207.204 <- rep(0,nt)
    sk.238.204 <- rep(0,nt)
    ti <- rep(0,nt)
    young <- which(tt < 3700)
    old <- which(tt >= 3700)
    sk.206.204[young] <- 11.152
    sk.207.204[young] <- 12.998
    sk.238.204[young] <- 9.74
    ti[young] <- 3700
    sk.206.204[old] <- 9.307
    sk.207.204[old] <- 10.294
    sk.238.204[old] <- 7.19
    ti[old] <- 4570
    U238U235 <- settings('iratio','U238U235')[1]
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    i64 <- sk.206.204 + sk.238.204*(exp(l8*ti)-exp(l8*tt))
    i74 <- sk.207.204 + sk.238.204*(exp(l5*ti)-exp(l5*tt))/U238U235
    if (inverse) out <- cbind(1/i64,i74/i64)
    else out <- cbind(i64,i74)
    out
}
