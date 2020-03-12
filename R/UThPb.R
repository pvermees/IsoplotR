ludwig.UThPb <- function(x,exterr=FALSE,alpha=0.05,model=1,
                         anchor=list(FALSE,NA)){
    lta0b0 <- get.lta0b0.init(x,anchor=anchor)
    fit <- get.lta0b0c0.UThPb(x,lta0b0=lta0b0,exterr=exterr,
                              model=model,anchor=anchor)
    out <- exponentiate_ludwig(fit,format=x$format)
    out$n <- length(x)
    out
}

get.lta0b0c0.UThPb <- function(x,lta0b0,exterr=FALSE,model=1,
                               anchor=list(FALSE,NA),...){
    lt <- lta0b0[1]
    a0 <- lta0b0[2]
    b0 <- lta0b0[3]
    NR <- 4
    ns <- length(x)
    zeros <- rep(0,ns)
    XYZW <- list(X=zeros,Y=zeros,Z=zeros,W=zeros,
                 E=matrix(0,NR*ns,NR*ns))
    for (i in 1:ns){
        wd <- wetherill(x,i=i)
        XYZW$X[i] <- wd$x['Pb207U235']
        XYZW$Y[i] <- wd$x['Pb206U238']
        XYZW$Z[i] <- wd$x['Pb208Th232']
        XYZW$W[i] <- wd$x['Th232U238']
        XYZW$E[(0:(NR-1))*ns+i,(0:(NR-1))*ns+i] <- wd$cov
    }
    c0 <- log(XYZW$Z - exp(lambda('Th232')[1]*exp(lt)) + 1)
    if (model==3){
        init <- c(lt,a0,b0,0,c0)
        names(init) <- c('log(t)','log(a0)','log(b0)',
                         'log(w)',paste0('log(c0[',1:ns,'])'))
        lower <- c(lt,a0,b0,-8,c0)-2
        upper <- c(lt,a0,b0,0,c0)+2
        NP <- 4
    } else {
        init <- c(lt,a0,b0,c0)
        names(init) <- c('log(t)','log(a0)','log(b0)',
                         paste0('log(c0[',1:ns,'])'))
        lower <- init-1
        upper <- init+2
        NP <- 3
    }
    fixed <- rep(FALSE,length(init))
    fit <- optifix(parms=init,fn=LL.lud.UThPb,XYZW=XYZW,
                   method="L-BFGS-B",x=x,exterr=exterr,
                   fixed=fixed,lower=lower,upper=upper,
                   hessian=TRUE,LL=TRUE,
                   control=list(fnscale=-1),...)
    out <- LL.lud.UThPb(lta0b0wc0=fit$par,x=x,LL=FALSE,XYZW=XYZW)
    out$model <- model
    out$logpar <- fit$par[1:NP]
    out$logcov <- solve(-fit$hessian)[1:NP,1:NP]
    out
}

LL.lud.UThPb <- function(lta0b0wc0,x,exterr=FALSE,LL=TRUE,XYZW){
    ns <- length(x)
    NP <- length(lta0b0wc0)
    tt <- exp(lta0b0wc0[1])
    a0 <- exp(lta0b0wc0[2])
    b0 <- exp(lta0b0wc0[3])
    model3 <- NP>(ns+3)
    if (model3){
        w <- exp(lta0b0wc0[4])
        c0 <- exp(lta0b0wc0[4+(1:ns)])
    } else {
        c0 <- exp(lta0b0wc0[3+(1:ns)])
    }
    X <- XYZW$X
    Y <- XYZW$Y
    Z <- XYZW$Z
    W <- XYZW$W
    E <- XYZW$E
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
        ED <- ED + Jw%*%w^2%*% t(Jw)
    }
    D <- c(K,L,M)
    SS <- D %*% solve(ED) %*% D
    if (LL){
        detED <- determinant(ED,logarithm=TRUE)$modulus
        out <- -(NP*ns*log(2*pi) + detED + SS)/2
    } else {
        out <- list()
        out$df <- 3*ns-NP
        out$mswd <- as.vector(SS/out$df)
        out$p.value <- as.numeric(1-stats::pchisq(SS,out$df))
    }
    out
}

init.lw.UThPb <- function(lt,a0,b0,c0,x,XYZW){
    LL <- function(lw,lt,a0,b0,c0,x,XYZW){
        LL.lud.UThPb(lta0b0wc0=c(lt,a0,b0,lw,c0),
                     x=x,XYZW=XYZW,LL=TRUE)
    }
    stats::optimise(LL,interval=c(-10,2),lt=lt,a0=a0,b0=b0,
                    c0=c0,x=x,XYZW=XYZW,maximum=TRUE)$maximum
}
