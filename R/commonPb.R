# option = 1: Stacey-Kramers
# option = 2: isochron
# option = 3: assumed Pb-composition
common.Pb.correction <- function(x,...){ UseMethod("common.Pb.correction",x) }
common.Pb.correction.default <- function(x,method='common.Pb.correction',...){
    stop('No default method for common.Pb.correction function')
}
common.Pb.correction.UPb <- function(x,option=1){
    ns <- length(x)
    if (option == 1) out <- common.Pb.stacey.kramers.UPb(x)
    else if (option == 2) out <- common.Pb.isochron.UPb(x)
    else if (option == 3) out <- common.Pb.nominal.UPb(x)
    else out <- x
    out$x.raw <- x$x
    out
}
common.Pb.correction.PbPb <- function(x,option=1){
    ns <- length(x)
    if (option == 1) out <- common.Pb.stacey.kramers.PbPb(x)
    else if (option == 2) out <- common.Pb.isochron.PbPb(x)
    else if (option == 3) out <- common.Pb.nominal.PbPb(x)
    else out <- x
    out$x.raw <- x$x
    out
}

common.Pb.isochron.UPb <- function(x){
    if (x$format<4){
        d <- data2york(x,wetherill=FALSE)
        fit <- york(d)
        i76 <- fit$a[1]
    } else {
        lud <- ludwig(x)
        i64 <- lud$par['64i']
        i74 <- lud$par['74i']
    }
    if (x$format < 4)
        out <- Pb.correction.without.204(x,i76)
    else
        out <- Pb.correction.with.204(x,i64,i74)
    out
}
common.Pb.isochron.PbPb <- function(x){
    out <- x
    fit <- isochron(x,inverse=FALSE,plot=FALSE)
    i74 <- fit$a[1]
    i64 <- fit$a[1]/fit$b[1]
    Pb.correction.for.PbPb(x,i64,i74)
}

common.Pb.stacey.kramers.UPb <- function(x){
    tt <- 1000
    for (i in 1:5){
        i6474 <- stacey.kramers(tt)
        i64 <- i6474[1]
        i74 <- i6474[2]
        if (x$format < 4)
            out <- Pb.correction.without.204(x,i74/i64)
        else
            out <- Pb.correction.with.204(x,i64,i74)
        tt <- get.Pb206U238.age(out)[,1]
    }
    out
}
common.Pb.stacey.kramers.PbPb <- function(x){
    tt <- 1000
    for (i in 1:5){
        i6474 <- stacey.kramers(tt)
        i64 <- i6474[1]
        i74 <- i6474[2]
        out <- Pb.correction.for.PbPb(x,i64,i74)
        tt <- PbPb.age(out)[,1]        
    }
    out
}

common.Pb.nominal.UPb <- function(x){
    i64 <- settings('iratio','Pb206Pb204')[1]
    i74 <- settings('iratio','Pb207Pb204')[1]
    if (x$format < 4)
        out <- Pb.correction.without.204(x,i74/i64)
    else
        out <- Pb.correction.with.204(x,i64,i74)
    out
}
common.Pb.nominal.PbPb <- function(x){
    out <- x
    i64 <- settings('iratio','Pb206Pb204')[1]
    i74 <- settings('iratio','Pb207Pb204')[1]
    Pb.correction.for.PbPb(x,i74,i64)
}

Pb.correction.without.204 <- function(x,i76){
    ns <- length(x)
    ni <- length(i76)
    out <- x
    m76 <- get.Pb207Pb206.ratios(x)[,1]
    m86 <- get.U238Pb206.ratios(x)[,1]
    tint <- rep(0,ns) # intercept age
    for (i in 1:ns){
        if (ni>1)
            tint[i] <- project.concordia(m76[i],m86[i],i76[i])
        else
            tint[i] <- project.concordia(m76[i],m86[i],i76)
    }
    if (x$format == 1){
        out$x[,'Pb207U235'] <- age_to_Pb207U235_ratio(tint)[,'75']
        out$x[,'errPb207U235'] <- x$x[,'errPb207U235']
        out$x[,'Pb206U238'] <- age_to_Pb206U238_ratio(tint)[,'68']
        out$x[,'errPb206U238'] <- x$x[,'errPb206U238']
        out$x[,'rhoXY'] <- x$x[,'rhoXY']        
    } else if (x$format == 2){
        out$x[,'U238Pb206'] <- age_to_U238Pb206_ratio(tint)[,'86']
        out$x[,'errU238Pb206'] <- x$x[,'errU238Pb206']
        out$x[,'Pb207Pb206'] <- age_to_Pb207Pb206_ratio(tint)[,'76']
        out$x[,'errPb207Pb206'] <- x$x[,'errPb207Pb206']
        out$x[,'rhoXY'] <- x$x[,'rhoXY']
    } else if (x$format == 3){
        out$x[,'Pb207Pb206'] <- age_to_Pb207Pb206_ratio(tint)[,'76']
        out$x[,'errPb207Pb206'] <- x$x[,'errPb207Pb206']
        out$x[,'Pb207U235'] <- age_to_Pb207U235_ratio(tint)[,'75']
        out$x[,'errPb207U235'] <- x$x[,'errPb207U235']
        out$x[,'Pb206U238'] <- age_to_Pb206U238_ratio(tint)[,'68']
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
        colnames(out$x) <- c('Pb207U235','errPb207U235','Pb206U238','errPb206U238','rhoXY')
        out$x[,'Pb207U235'] <- x$x[,'Pb207U235'] - i74*x$x[,'Pb204U238']*U238U235
        out$x[,'errPb207U235'] <- x$x[,'errPb207U235']
        out$x[,'Pb206U238'] <- x$x[,'Pb206U238'] - i64*x$x[,'Pb204U238']
        out$x[,'errPb206U238'] <- x$x[,'errPb206U238']
        out$x[,'rhoXY'] <- x$x[,'rhoXY']
        out$format <- 1
    } else if (x$format == 5){
        out$x <- matrix(0,ns,5)
        colnames(out$x) <- c('U238Pb206','errU238Pb206','Pb207Pb206','errPb207Pb206','rhoXY')
        out$x[,'U238Pb206'] <- x$x[,'U238Pb206']/(1 - i64*x$x[,'Pb204Pb206'])
        out$x[,'errU238Pb206'] <- x$x[,'errU238Pb206']
        out$x[,'Pb207Pb206'] <-
            (x$x[,'Pb207Pb206'] - x$x[,'Pb204Pb206']*i74)/(1 - x$x[,'Pb204Pb206']*i64)
        out$x[,'errPb207Pb206'] <- x$x[,'errPb207Pb206']
        out$x[,'rhoXY'] <- x$x[,'rhoXY']
        out$format <- 2
    } else if (x$format == 6){
        out$x <- matrix(0,ns,9)
        colnames(out$x) <- c('Pb207Pb206','errPb207Pb206',
                             'Pb207U235','errPb207U235',
                             'Pb206U238','errPb206U238',
                             'rhoXY','rhoXZ','rhoYZ')
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
Pb.correction.for.PbPb <- function(x,i64,i74){
    out <- x
    if (x$format %in% c(1,3)){
        m64 <- x$x[,'Pb206Pb204']
        m74 <- x$x[,'Pb207Pb204']
        r64 <- m64 - i64
        r74 <- m74 - i74
        out$x[,'Pb206Pb204'] <- r64
        out$x[,'Pb207Pb204'] <- r74
    }
    if (x$format == 2){
        m46 <- x$x[,'Pb204Pb206']
        m76 <- x$x[,'Pb207Pb206']
        r64 <- 1/m46 - i64
        r46 <- 1/r64
        r74 <- m76/m46 - i74
        r76 <- r74*r46
        out$x[,'Pb204Pb206'] <- r46
        out$x[,'Pb207Pb206'] <- r76
        if (FALSE){
            dr46.dr64 <- -1/r64^2
            dr64.dm46 <- -1/m46^2
            dr74.dm46 <- -m76/m46^2
            dr74.dm76 <- 1/m46
            dr46.dm46 <- dr46.dr64*dr64.dm46
            dr46.dm76 <- 0
            J11 <- dr46.dm46                     # dr46.dm46
            J12 <- dr46.dm76                     # dr46.dm76
            J21 <- dr74.dm46*r46 + r74*dr46.dm46 # dr76.dm46
            J22 <- dr74.dm76*r46 + r74*dr46.dm76 # dr76.dm76
            E11 <- x$x[,'errPb204Pb206']^2
            E22 <- x$x[,'errPb207Pb206']^2
            E12 <- x$x[,'rho']*
                x$x[,'errPb204Pb206']*x$x[,'errPb207Pb206']
            covmat <- errorprop(J11,J12,J21,J22,E11,E22,E12)
            out$x[,'errPb204Pb206'] <- sqrt(covmat[,'varX'])
            out$x[,'errPb207Pb206'] <- sqrt(covmat[,'varY'])
            out$x[,'rho'] <- covmat[,'cov']/
                sqrt(covmat[,'varX']*covmat[,'varY'])
        }
    } else if (x$format == 3){
        dat <- PbPb.normal.ratios(x)
        r76 <- r74/r64
        dr74.dm74 <- 1
        dr64.dm74 <- 0
        dr74.dm64 <- 0
        dr64.dm64 <- 1
        dr76.dm64 <- (r64*dr74.dm64 - r74*dr64.dm64)/r64^2
        dr76.dm74 <- (r64*dr74.dm74 - r74*dr64.dm74)/r64^2
        if (FALSE){
            J <- matrix(0,1,2)
            E <- matrix(0,2,2)
            for (i in 1:length(x)){
                J[1,1] <- dr76.dm64[i]
                J[1,2] <- dr76.dm74[i]
                E[1,1] <- dat[i,'errPb206Pb204']^2
                E[2,2] <- dat[i,'errPb207Pb204']^2
                E[1,2] <- dat[i,'rho']*dat[i,'errPb206Pb204']*dat[i,'errPb207Pb204']
                E[2,1] <- E[1,2]
                out$x[i,'errPb207Pb206'] <- sqrt(J %*% E %*% t(J))
            }
        }
        out$x[,'Pb207Pb206'] <- r76
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
