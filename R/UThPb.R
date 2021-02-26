ludwig.UThPb <- function(x,exterr=FALSE,alpha=0.05,model=1,anchor=0){
    init <- get.lta0b0wc0.init(x=x,anchor=anchor,model=model)
    fit <- get.lta0b0wc0.UThPb(x,init=init,exterr=exterr,
                               model=model,anchor=anchor)
    out <- exponentiate_ludwig(fit,format=x$format)
    out$n <- length(x)
    out
}

get.lta0b0wc0.init <- function(x,model=1,anchor=0){
    ns <- length(x)
    if (model==3){
        lta0b0w <- get.lta0b0w(x=x,model=3,anchor=anchor)$logpar
        c0 <- log(data2ludwig(x=x,lta0b0w=lta0b0w)$c0)
        out <- c(lta0b0w,c0)
        names(out) <- c('log(t)','log(a0)','log(b0)',
                        'log(w)',paste0('log(c0[',1:ns,'])'))
    } else {
        lta0b0 <- get.lta0b0.init(x,anchor=anchor)
        c0 <- log(data2ludwig(x=x,lta0b0w=lta0b0)$c0)
        out <- c(lta0b0,c0)
        names(out) <- c('log(t)','log(a0)','log(b0)',
                        paste0('log(c0[',1:ns,'])'))
    }
    out
}

get.lta0b0wc0.UThPb <- function(x,init,exterr=FALSE,model=1,anchor=0,...){
    XYZW <- get_XYZW(x)
    fixed <- rep(FALSE,length(init))
    lower <- init-2
    upper <- init+2
    if (model==3){
        lower[4] <- lower[1]-4
        upper[4] <- upper[1]
    }
    fit <- optifix(parms=init,fn=LL.lud.UThPb,XYZW=XYZW,
                   method="L-BFGS-B",lower=lower,upper=upper,
                   x=x,exterr=exterr,fixed=fixed,hessian=TRUE,
                   LL=TRUE,control=list(fnscale=-1),...)
    out <- LL.lud.UThPb(lta0b0wc0=fit$par,x=x,LL=FALSE,XYZW=XYZW)
    out$model <- model
    if (model==3) NP <- 4 else NP <- 3
    out$logpar <- fit$par[1:NP]
    out$logcov <- solve(-fit$hessian)[1:NP,1:NP]
    out
}

get_XYZW <- function(x){
    ns <- length(x)
    zeros <- rep(0,ns)
    out <- list(X=zeros,Y=zeros,Z=zeros,W=zeros,
                E=matrix(0,4*ns,4*ns))
    for (i in 1:ns){
        wd <- wetherill(x,i=i)
        out$X[i] <- wd$x['Pb207U235']
        out$Y[i] <- wd$x['Pb206U238']
        out$Z[i] <- wd$x['Pb208Th232']
        out$W[i] <- wd$x['Th232U238']
        out$E[(0:3)*ns+i,(0:3)*ns+i] <- wd$cov
    }
    out
}

LL.lud.UThPb <- function(lta0b0wc0,x,exterr=FALSE,LL=TRUE,XYZW){
    tt <- exp(lta0b0wc0[1])
    a0 <- exp(lta0b0wc0[2])
    b0 <- exp(lta0b0wc0[3])
    X <- XYZW$X
    Y <- XYZW$Y
    Z <- XYZW$Z
    W <- XYZW$W
    E <- XYZW$E
    ns <- length(x)
    NP <- length(lta0b0wc0)
    model3 <- NP>(ns+3)
    if (model3){
        w <- exp(lta0b0wc0[4])
        c0 <- exp(lta0b0wc0[4+(1:ns)])
    } else {
        c0 <- exp(lta0b0wc0[3+(1:ns)])
    }
    U <- iratio('U238U235')[1]
    D <- mclean(tt=tt,d=x$d,exterr=exterr)
    K <- X - b0*U*W*c0 - D$Pb207U235
    L <- Y - a0*W*c0 - D$Pb206U238
    M <- Z - c0 - D$Pb208Th232
    ns <- length(X)
    J <- matrix(0,3*ns,4*ns)
    i1 <- 1:ns
    i2 <- (ns+1):(2*ns)
    i3 <- (2*ns+1):(3*ns)
    i4 <- (3*ns+1):(4*ns)
    J[i1,i1] <- diag(ns)
    J[i2,i2] <- diag(ns)
    J[i3,i3] <- diag(ns)
    J[i1,i4] <- -U*b0*c0
    J[i2,i4] <- -a0*c0
    ED <- J %*% E %*% t(J)
    if (model3){
        Jw <- matrix(0,3*ns,1)
        Jw[i1] <- -D$dPb207U235dt
        Jw[i2] <- -D$dPb206U238dt
        Jw[i3] <- -D$dPb208Th232dt
        ED <- ED + Jw %*% (w^2) %*% t(Jw)
    }
    O <- blockinverse3x3(ED[i1,i1],ED[i1,i2],ED[i1,i3],
                         ED[i2,i1],ED[i2,i2],ED[i2,i3],
                         ED[i3,i1],ED[i3,i2],ED[i3,i3])
    D <- c(K,L,M)
    SS <- D %*% O %*% D
    if (LL){
        detED <- determinant(ED,logarithm=TRUE)$modulus
        out <- -(3*ns*log(2*pi) + detED + SS)/2
    } else {
        out <- list()
        out$df <- 3*ns-NP
        out$mswd <- as.vector(SS/out$df)
        out$p.value <- as.numeric(1-stats::pchisq(SS,out$df))
    }
    out
}
