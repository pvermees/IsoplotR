ludwig.UThPb <- function(x,exterr=FALSE,alpha=0.05,model=1,
                         anchor=list(FALSE,NA)){
    init <- get.lta0b0.init(x,anchor=anchor)
    init <- c(2.8710214,-0.9425673,-2.4630234)
    fit <- get.lta0b0c0w(x,lta0b0=init,exterr=exterr,
                         model=model,anchor=anchor)
    out <- exponentiate_ludwig(fit,format=x$format)
    out$n <- length(x)
    out
}

get.lta0b0c0w <- function(x,lta0b0,exterr=FALSE,model=1,
                          anchor=list(FALSE,NA),w=NA,...){
    NP <- length(lta0b0)
    lt <- lta0b0[1]
    a0 <- lta0b0[2]
    b0 <- lta0b0[3]
    NR <- 4
    ns <- length(x)
    zeros <- rep(0,ns)
    XYZW <- list(X=zeros,Y=zeros,Z=zeros,W=zeros,E=matrix(0,NR*ns,NR*ns))
    for (i in 1:ns){
        wd <- wetherill(x,i=i)
        XYZW$X[i] <- wd$x['Pb207U235']
        XYZW$Y[i] <- wd$x['Pb206U238']
        XYZW$Z[i] <- wd$x['Pb208Th232']
        XYZW$W[i] <- wd$x['Th232U238']
        XYZW$E[(0:(NR-1))*ns+i,(0:(NR-1))*ns+i] <- wd$cov
    }
    l <- data2ludwig(x,lta0b0)
    c0 <- l$c0
    init <- c(lt,a0,b0,c0)
    lower <- c(c(lt,a0,b0)-1,c0/2)
    upper <- c(c(lt,a0,b0)+2,c0*2)
    fixed <- rep(FALSE,length(init))
    fit <- optifix(parms=init,fn=LL.lud.UThPb,XYZW=XYZW,
                   method="L-BFGS-B",x=x,exterr=exterr,
                   fixed=fixed,lower=lower,upper=upper,
                   hessian=TRUE,...)
    out <- list()
    out$model <- model
    out$logpar <- fit$par[1:3]
    out$logcov <- solve(fit$hessian)[1:3,1:3]
    out$df <- 2*ns-3
    out$mswd <- fit$value/out$df
    out$p.value <- as.numeric(1-stats::pchisq(fit$value,out$df))
    out
}

LL.lud.UThPb <- function(lta0b0c0,x,exterr=FALSE,LL=TRUE,XYZW){
    tt <- exp(lta0b0c0[1])
    a0 <- exp(lta0b0c0[2])
    b0 <- exp(lta0b0c0[3])
    c0 <- lta0b0c0[-(1:3)]
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
    D <- c(K,L,M)
    SS <- D %*% solve(ED) %*% D
    SS
}

mswd.lud.UThPb <- function(lta0b0c0,x,anchor=list(FALSE,NA)){
    ns <- length(x)
    out <- list()
    LL.lud.UThPb(lta0b0c0,x,XYZW)
    if (out$df>0){
        out$mswd <- as.vector(SS/out$df)
        out$p.value <- as.numeric(1-stats::pchisq(SS,out$df))
    } else {
        out$mswd <- 1
        out$p.value <- 1
    }
    out
}
