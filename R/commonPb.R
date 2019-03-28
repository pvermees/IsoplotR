# option = 1: Stacey-Kramers
# option = 2: isochron
# option = 3: assumed Pb-composition
common.Pb.correction <- function(x,option=1,calcit=rep(TRUE,length(x))){
    ns <- length(x)
    out <- x
    out$x.raw <- out$x
    out$x <- matrix(0,ns,5)
    if (option == 1)
        out$x <- common.Pb.stacey.kramers(x)
    else if (option == 2)
        out$x <- common.Pb.isochron(x,calcit=calcit)
    else if (option == 3)
        out$x <- common.Pb.nominal(x)
    else out$x <- x
    if (x$format<4){
        out$format <- 2
        colnames(out$x) <- c('U238Pb206','errU238Pb206',
                             'Pb207Pb206','errPb207Pb206','rhoXY')
    } else {
        out$format <- 1
        colnames(out$x) <- c('Pb207U235','errPb207U235',
                             'Pb206U238','errPb206U238','rhoXY')
    }
    out
}

correct.common.Pb.without.204 <- function(x,i,c76,lower=TRUE){
    tw <- tera.wasserburg(x,i)
    m86 <- tw$x['U238Pb206']
    m76 <- tw$x['Pb207Pb206']
    tint <- project.concordia(m86,m76,c76,d=x$d,lower=lower)
    cctw <- age_to_terawasserburg_ratios(tt=tint,st=0,d=x$d)
    r86 <- cctw$x['U238Pb206']
    r76 <- cctw$x['Pb207Pb206']
    f <- (m76-r76)/(m76-c76)
    E <- tw$cov/((1-f)^2)
    sr86 <- sqrt(E[1,1])
    sr76 <- sqrt(E[2,2])
    rho <- stats::cov2cor(E)[1,2]
    out <- c(r86,sr86,r76,sr76,rho)
    out
}
correct.common.Pb.with.204 <- function(x,i,c46,c47){
    ir <- get.UPb.isochron.ratios(x,i) # 68, 46, 75, 47
    m68 <- ir$x['Pb206U238']
    m46 <- ir$x['Pb204Pb206']
    m75 <- ir$x['Pb207U235']
    m47 <- ir$x['Pb204Pb207']    
    r75 <- m75*(1-m47/c47)
    r68 <- m68*(1-m46/c46)
    J <- matrix(0,2,4)
    J[1,3] <- 1-m47/c47
    J[1,4] <- -m75/c47
    J[2,1] <- 1-m46/c46
    J[2,2] <- -m68/c46
    E <- J %*% ir$cov %*% t(J)
    sr75 <- sqrt(E[1,1])
    sr68 <- sqrt(E[2,2])
    rho <- stats::cov2cor(E)[1,2]
    out <- c(r75,sr75,r68,sr68,rho)
    out
}

common.Pb.stacey.kramers <- function(x){
    ns <- length(x)
    out <- matrix(0,ns,5)
    if (x$format < 4){
        for (i in 1:ns){
            tint <- stats::optimise(SS.SK.without.204,
                                    interval=c(0,5000),x=x,i=i)$minimum
            i6474 <- stacey.kramers(tint)
            c76 <- i6474[2]/i6474[1]
            out[i,] <- correct.common.Pb.without.204(x,i,c76,lower=FALSE)
        }
    } else {
        for (i in 1:ns){
            tint <- stats::optimise(SS.SK.with.204,
                                    interval=c(0,5000),x=x,i=i)$minimum
            c6474 <- stacey.kramers(tint)
            c46 <- 1/c6474[1]
            c47 <- 1/c6474[2]
            out[i,] <- correct.common.Pb.with.204(x,i,c46,c47)
        }
    }
    out
}

common.Pb.isochron <- function(x,calcit=rep(TRUE,length(x))){
    fit <- ludwig(subset(x,subset=calcit))
    ns <- length(x)
    out <- matrix(0,ns,5)
    tt <- fit$par[1]
    if (x$format<4){
        rr <- age_to_terawasserburg_ratios(tt,d=x$d)$x
        slope <- (rr['Pb207Pb206']-fit$par['76i'])/rr['U238Pb206']
        m76 <- get.Pb207Pb206.ratios(x)[,1]
        m86 <- get.U238Pb206.ratios(x)[,1]
        c76 <- m76 - slope*m86
        for (i in 1:ns){
            out[i,] <- correct.common.Pb.without.204(x,i,c76[i],lower=TRUE)
        }
    } else {
        r68 <- age_to_Pb206U238_ratio(tt=tt,st=0,d=x$d)
        slope.68 <- -r68[1]/fit$par['64i']
        r75 <- age_to_Pb207U235_ratio(tt=tt,st=0,d=x$d)
        slope.75 <- -r75[1]/fit$par['74i']
        for (i in 1:ns){
            rr <- get.UPb.isochron.ratios(x,i)
            m46 <- rr$x['Pb204Pb206']
            m68 <- rr$x['Pb206U238']
            c46 <- m46 - slope.68/m68
            m47 <- rr$x['Pb204Pb207']
            m75 <- rr$x['Pb207U235']
            c47 <- m47 - slope.75/m75
            out[i,] <- correct.common.Pb.with.204(x,i,c46,c47)
        }
    }
    out
}

common.Pb.nominal <- function(x){
    ns <- length(x)
    out <- matrix(0,ns,5)
    if (x$format < 4){
        c76 <- settings('iratio','Pb207Pb206')[1]
        for (i in 1:ns){
            out[i,] <- correct.common.Pb.without.204(x,i,c76,lower=TRUE)
        }
    } else {
        c46 <- 1/settings('iratio','Pb206Pb204')[1]
        c47 <- 1/settings('iratio','Pb207Pb204')[1]
        for (i in 1:ns){
            out[i,] <- correct.common.Pb.with.204(x,i,c46,c47)
        }
    }
    out
}

SS.SK.without.204 <- function(tt,x,i){
    tw <- tera.wasserburg(x,i)
    X <- tw$x['U238Pb206']
    Y <- tw$x['Pb207Pb206']
    i6474 <- stacey.kramers(tt)
    cct <- age_to_terawasserburg_ratios(tt,st=0,d=x$d)
    a <- i6474[2]/i6474[1] # intercept
    b <- (cct$x['Pb207Pb206']-a)/cct$x['U238Pb206'] # slope
    omega <- solve(tw$cov)
    x.fitted <- (X*omega[1,1]+Y*omega[1,2]-a*omega[1,2])/(omega[1,1]+b*omega[1,2])
    y.fitted <- a + b*x.fitted
    d <- cbind(X-x.fitted,Y-y.fitted)
    as.numeric(d %*% omega %*% t(d))
}
SS.SK.with.204 <- function(tt,x,i){
    wi <- wetherill(x,i=i)
    i6474 <- stacey.kramers(tt)
    i64 <- i6474[1]
    i74 <- i6474[2]
    U <- iratio('U238U235')[1]
    ccw <- list(x=rep(0,2),cov=matrix(0,2,2))
    ccw$x[1] <- wi$x['Pb207U235'] - i74*wi$x['Pb204U238']*U
    ccw$x[2] <- wi$x['Pb206U238'] - i64*wi$x['Pb204U238']
    J <- matrix(0,2,3)
    J[1,1] <- 1
    J[1,3] <- -i74*U
    J[2,2] <- 1
    J[2,3] <- -i64
    ccw$cov <- J %*% wi$cov %*% t(J)
    LL.concordia.age(tt,ccw,mswd=FALSE,exterr=FALSE,d=x$d)
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
