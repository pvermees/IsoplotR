ludwig2d <- function(x,type=1,model=1,anchor=0,exterr=FALSE){
    if (model==2){
        fit <- ludwig2d_model2(x=x,type=type,anchor=anchor,exterr=exterr)
    } else {
        fit <- ludwig2d_helper(x=x,model=model,type=type,
                               anchor=anchor,exterr=exterr)
    }

    out <- fit
    np <- ifelse(model==3,4,3)
    out$par <- rep(NA,np)
    out$cov <- matrix(0,np,np)
    out$logpar <- rep(NA,np)
    out$logcov <- matrix(0,np,np)
    
    if (type%in%c(1,3)) i <- c(1,2)
    else if (type%in%c(2,4)) i <- c(1,3)
    else stop('Invalid isochron type.')
    if (model==3) i <- c(i,4)

    out$par[i] <- fit$par
    out$cov[i,i] <- fit$cov
    out$logpar[i] <- fit$logpar
    out$logcov[i,i] <- fit$logcov
    
    pnames <- c('t','a0','b0')
    if (model==3) pnames <- c(pnames,'w')
    names(out$par) <- pnames
    rownames(out$cov) <- pnames
    colnames(out$cov) <- pnames

    lnames <- c('log(t)','log(a0)','log(b0)')
    if (model==3) lnames <- c(lnames,'log(w)')
    names(out$logpar) <- lnames
    rownames(out$logcov) <- lnames
    colnames(out$logcov) <- lnames

    out$model <- model
    out$n <- length(x)
    out
}

ludwig2d_model2 <- function(x,type=1,anchor=0,exterr=FALSE){

    out <- list(par=rep(NA,2),cov=matrix(0,2,2),
                logpar=rep(NA,2),logcov=matrix(0,2,2))
        
    if (x$format%in%(4:6)){
        
        XY <- data2york(x,option=(3:4)[type])
        if (anchor[1]<1){
            fit <- stats::lm( XY[,'Y'] ~ XY[,'X'] )
            ab <- fit$coefficients
            covmat <- stats::vcov(fit)
        } else if (anchor[1]==1) {
            y0 <- 1/ifelse(type==1,iratio('Pb206Pb204')[1],iratio('Pb207Pb204')[1])
            fit <- stats::lm( I(XY[,'Y']-y0) ~ 0 + XY[,'X'] )
            ab <- c(y0,fit$coefficients)
            covmat <- matrix(0,2,2)
            covmat[2,2] <- stats::vcov(fit)
        } else {
            D <- mclean(anchor[2],d=x$d)
            x0 <- 1/ifelse(type==1,D$Pb206U238,D$Pb207U235)
            fit <- stats::lm( I(XY[,'Y']-0) ~ 0 + I(XY[,'X']-x0) )
            y0 <- -x0*fit$coefficients
            ab <- c(y0,fit$coefficients)
            covmat <- matrix(0,2,2)
            covmat[2,2] <- stats::vcov(fit)
        }
        DP <- -ab[2]/ab[1]
        Dd <- 1/ab[1]
        if (type==1){
            tt <- get.Pb206U238.age(DP)[1]
            D <- mclean(tt,d=x$d,exterr=exterr)
            dtdDP <- 1/D$dPb206U238dt
            E <- matrix(0,6,6)
            E[3:6,3:6] <- getEl('U238')
            J <- matrix(0,2,6)
            J[1,3] <- dtdDP*D$dPb206U238dl38
            J[1,4] <- dtdDP*D$dPb206U238dl34
            J[1,5] <- dtdDP*D$dPb206U238dl30
            J[1,6] <- dtdDP*D$dPb206U238dl26
        } else {
            tt <- get.Pb207U235.age(DP)[1]
            D <- mclean(tt,d=x$d,exterr=exterr)
            dtdDP <- 1/D$dPb207U235dt
            E <- matrix(0,4,4)
            E[3:4,3:4] <- getEl('U235')
            J <- matrix(0,2,4)
            J[1,3] <- dtdDP*D$dPb207U235dl35
            J[1,4] <- dtdDP*D$dPb207U235dl31
        }
        E[1:2,1:2] <- covmat
        J[1,1] <- dtdDP*ab[2]/ab[1]^2
        J[1,2] <- -dtdDP/ab[1]
        J[2,1] <- -1/ab[1]^2
        out$par <- c(tt,Dd)
        out$cov <- J%*%E%*%t(J)
        out$logpar <- log(out$par)
        J <- diag(1/out$par)
        out$logcov <- J%*%out$cov%*%t(J)
        
    } else if (x$format%in%(7:8)){
        
        LL <- function(lta0,x,np=2,option=1,exterr=FALSE){
            ns <- length(x)
            tt <- exp(lta0[1])
            a0 <- exp(lta0[2])
            df <- ns-np
            D <- mclean(tt,d=x$d,exterr=exterr)
            if (option==6){
                x0 <- 1/D$Pb206U238
                y0 <- 1/a0
            } else if (option==7){
                x0 <- 1/D$Pb207U235
                y0 <- 1/a0
            } else {
                x0 <- 1/D$Pb208Th232
                y0 <- a0
            }
            XY <- data2york(x,option=option,tt=tt)
            yp <- y0*(1-XY[,'X']/x0)
            SS <- sum((yp-XY[,'Y'])^2)
            if (exterr){
                dypdx0 <- y0*XY[,'X',drop=FALSE]/x0^2
                dx0dDP <- -x0^2
                dDPdl <- matrix(0,1,7)
                if (option==6){
                    dDPdl[1] <- D$dPb206U238dl38
                    dDPdl[3] <- D$dPb206U238dl34
                    dDPdl[6] <- D$dPb206U238dl30
                    dDPdl[7] <- D$dPb206U238dl26
                } else if (option==7){
                    dDPdl[2] <- D$dPb207U235dl35
                    dDPdl[5] <- D$dPb207U235dl31
                } else {
                    dDPdl[4] <- D$dPb208Th232dl32
                }
                J <- (dypdx0 * dx0dDP) %*% dDPdl
                covmat <- diag(SS/df,ns,ns) + J%*%getEl()%*%t(J)
                out <- LL.norm(yp-XY[,'Y'],covmat)
            } else {
                out <- SS2LL(SS,ns,df)
            } 
            out
        }
        
        model2init <- function(x,option=1,anchor=0){
            init <- c(0,0)
            if (anchor[1]==2){
                ti <- anchor[2]
            } else {
                ti <- min(get.Pb206U238.age(x)[,1]) # first stab
                XY <- data2york(x,option=option,tt=ti)
                if (option==6){
                    init[1] <- min(get.Pb206U238.age(1/XY[,'X'])[,1])
                } else if (option==7){
                    init[1] <- min(get.Pb207U235.age(1/XY[,'X'])[,1])
                } else {
                    init[1] <- min(get.Pb208Th232.age(1/XY[,'X'])[,1])
                }
            }
            init[1] <- ti
            if (anchor[1]==1) {
                if (option==6) {
                    init[2] <- 1/iratio('Pb208Pb206')[1]
                } else if (option==7) {
                    init[2] <- 1/iratio('Pb208Pb207')[1]
                } else if (option==7) {
                    init[2] <- iratio('Pb208Pb206')[1]
                } else if (option==8) {
                    init[2] <- iratio('Pb208Pb207')[1]
                } else {
                    stop("Illegal isochron option")
                }
            } else {
                XY <- data2york(x,option=option,tt=0)
                init[2] <- stats::median(XY[,'Y'])
            }
            log(init)
        }
        
        option <- (6:9)[type]
        init <- model2init(x,option=option,anchor=anchor)
        fixed <- fixit(x=x,anchor=anchor,model=2,joint=FALSE)
        lower <- (init-5)[!fixed]
        upper <- (init+5)[!fixed]
        fit <- optifix(init,fixed=fixed,fn=LL,x=x,option=option,
                       method='L-BFGS-B',lower=lower,upper=upper,
                       exterr=exterr,hessian=TRUE)
        out$logpar <- fit$par
        out$logcov[!fixed,!fixed] <- solve(fit$hessian)
        out$par <- exp(out$logpar)
        J <- diag(out$par)
        out$cov <- J %*% out$logcov %*% t(J)
        
    } else {
        stop('2D ludwig regression is not available for this format')
    }
    out
} # end of ludwig2d_model2

ludwig2d_helper <- function(x,model=1,type=1,anchor=0,exterr=FALSE){

    fixed <- fixit(x=x,anchor=anchor,model=model,joint=FALSE)
    np <- length(fixed)
    out <- list(par=rep(NA,np),cov=matrix(0,np,np),
                logpar=rep(NA,np),logcov=matrix(0,np,np))
    
    if (model==1){
        model2fit <- ludwig2d_model2(x=x,type=type,anchor=anchor)
        init <- model2fit$logpar
    } else {
        model1fit <- ludwig2d_helper(x=x,type=type,anchor=anchor,exterr=exterr)
        LL0 <- LL_lud2d(lta0w=c(model1fit$logpar,-5),x=x,type=type)
        LLw <- LL_lud2d(lta0w=c(model1fit$logpar,-4),x=x,type=type)
        if (LLw>LL0){ # zero dispersion
            out$par[1:2] <- model1fit$par
            out$cov[1:2,1:2] <- model1fit$cov
            out$logpar <- c(model1fit$logpar,-Inf)
            out$logcov[1:2,1:2] <- model1fit$logcov
            return(out)
        } else {
            winit <- initwlud2d(x=x,lta0=model1fit$logpar,type=type)
            init <- c(model1fit$logpar,winit)
        }
    }
    lower <- (init-2)[!fixed]
    upper <- (init+2)[!fixed]
    fit <- optifix(init,fixed=fixed,LL_lud2d,method='L-BFGS-B',
                   lower=lower,upper=upper,exterr=exterr,
                   x=x,type=type,hessian=TRUE)
    out$logpar <- fit$par
    out$logcov[!fixed,!fixed] <- solve(fit$hessian)
    out$par <- exp(out$logpar)
    out$cov <- diag(out$par)%*%out$logcov%*%diag(out$par)
    
    out$df <- ifelse(anchor[1]<1,length(x)-2,length(x)-1)
    SS <- LL_lud2d(fit$par[1:2],x=x,type=type,exterr=FALSE,LL=FALSE)
    out$mswd <- SS/out$df
    out$p.value <- as.numeric(1-stats::pchisq(SS,out$df))
    out
} # end of ludwig2d_helper

initwlud2d <- function(x=x,lta0,type=1){
    LL <- function(w,lta0,x,type=1){
        LL_lud2d(c(lta0,w),x=x,type=type)
    }
    stats::optimise(LL,lta0=lta0,x=x,type=type,lower=-5,upper=5)$min
}

LL_lud2d <- function(lta0w,x,type=1,LL=TRUE,exterr=FALSE){
    np <- length(lta0w)
    tt <- exp(lta0w[1])
    a0 <- exp(lta0w[2])
    w <- ifelse(np>2,exp(lta0w[3]),0)
    ns <- length(x)
    Y <- rep(NA,ns)
    Z <- rep(NA,ns)
    A <- diag(1,ns,ns)
    U <- iratio('U238U235')[1]
    E <- matrix(0,2*ns+7,2*ns+7)
    J <- matrix(0,2*ns,2*ns+7)
    J[1:(2*ns),1:(2*ns)] <- diag(2*ns)
    D <- mclean(tt=tt,d=x$d,exterr=exterr)
    nc <- length(D$ThUi) # nc>1 if each aliquot has its own diseq correction
    j <- 1
    for (i in 1:ns){
        wd <- wetherill(x,i=i)
        if (nc>1) j <- i
        if (x$format<7){
            if (type==1){
                Y[i] <- wd$x['Pb206U238'] - D$Pb206U238
                Z[i] <- wd$x['Pb204U238']
                A[i,i] <- a0
                E[i,i] <- wd$cov['Pb206U238','Pb206U238']
                E[ns+i,ns+i] <- wd$cov['Pb204U238','Pb204U238']
                E[i,ns+i] <- wd$cov['Pb206U238','Pb204U238']
                E[ns+i,i] <- E[i,ns+i]
                J[i,2*ns+1] <- -D$dPb206U238dl38[j]     #dLdl38
                J[i,2*ns+3] <- -D$dPb206U238dl34[j]     #dLdl34
                J[i,2*ns+6] <- -D$dPb206U238dl30[j]     #dLdl30
                J[i,2*ns+7] <- -D$dPb206U238dl26[j]     #dLdl26
            } else if (type==2){
                Y[i] <- wd$x['Pb207U235'] - D$Pb207U235
                Z[i] <- wd$x['Pb204U238']*U
                A[i,i] <- a0
                E[i,i] <- wd$cov['Pb207U235','Pb207U235']
                E[ns+i,ns+i] <- wd$cov['Pb204U238','Pb204U238']*U^2
                E[i,ns+i] <- wd$cov['Pb207U235','Pb204U238']*U
                E[ns+i,i] <- E[i,ns+i]
                J[i,2*ns+2] <- -D$dPb207U235dl35[j]     #dKdl35
                J[i,2*ns+5] <- -D$dPb207U235dl31[j]     #dKdl31
            } else {
                stop("Invalid isochron type.")
            }
        } else {
            if (type%in%c(1,3)){
                Y[i] <- wd$x['Pb206U238'] - D$Pb206U238
                Z[i] <- wd$x['Pb208Th232'] - D$Pb208Th232
                A[i,i] <- a0*wd$x['Th232U238']
                E[i,i] <- wd$cov['Pb206U238','Pb206U238']
                E[ns+i,ns+i] <- wd$cov['Pb208Th232','Pb208Th232']
                E[i,ns+i] <- wd$cov['Pb206U238','Pb208Th232']
                E[ns+i,i] <- E[i,ns+i]
                J[i,2*ns+1] <- -D$dPb206U238dl38[j]     #dLdl38
                J[i,2*ns+3] <- -D$dPb206U238dl34[j]     #dLdl34
                J[i,2*ns+6] <- -D$dPb206U238dl30[j]     #dLdl30
                J[i,2*ns+7] <- -D$dPb206U238dl26[j]     #dLdl26
                J[ns+i,2*ns+4] <- -D$dPb208Th232dl32[j] #dMdl32
            } else if (type%in%c(2,4)){
                Y[i] <- wd$x['Pb207U235'] - D$Pb207U235
                Z[i] <- wd$x['Pb208Th232'] - D$Pb208Th232
                A[i,i] <- a0*wd$x['Th232U238']*U
                E[i,i] <- wd$cov['Pb207U235','Pb207U235']
                E[ns+i,ns+i] <- wd$cov['Pb208Th232','Pb208Th232']
                E[i,ns+i] <- wd$cov['Pb207U235','Pb208Th232']
                E[ns+i,i] <- E[i,ns+i]
                J[i,2*ns+2] <- -D$dPb207U235dl35[j]     #dKdl35
                J[i,2*ns+5] <- -D$dPb207U235dl31[j]     #dKdl31
                J[ns+i,2*ns+4] <- -D$dPb208Th232dl32[j] #dMdl32
            } else {
                stop("Invalid isochron type.")
            }
        }
    }
    E[(2*ns)+(1:7),(2*ns)+(1:7)] <- getEl()
    ED <- J%*%E%*%t(J)
    i1 <- 1:ns
    i2 <- ns+(1:ns)
    Jw <- matrix(0,2*ns,ns)
    diag(Jw[i1,i1]) <- -ifelse(type%in%c(1,3),D$dPb206U238dt,D$Pb207U235)
    if (type%in%c(2,4)) diag(Jw[i2,i1]) <- -D$Pb208Th232
    Ew <- Jw%*%t(Jw)*w^2 + ED
    O <- blockinverse(Ew[i1,i1],Ew[i1,i2],Ew[i1,i2],Ew[i2,i2],doall=TRUE)
    AA <- Y%*%O[i1,i1]%*%Y + Y%*%O[i1,i2]%*%Z + Z%*%O[i2,i1]%*%Y + Z%*%O[i2,i2]%*%Z
    BB <- Y%*%O[i1,i1]%*%A + Y%*%O[i1,i2] + Z%*%O[i2,i1]%*%A + Z%*%O[i2,i2]
    CC <- A%*%O[i1,i1]%*%Y + A%*%O[i1,i2]%*%Z + O[i2,i1]%*%Y + O[i2,i2]%*%Z
    DD <- A%*%O[i1,i1]%*%A + A%*%O[i1,i2] + O[i2,i1]%*%A + O[i2,i2]
    z <- as.vector(solve(DD+t(DD),t(BB)+CC))
    SS <- AA - BB%*%z - t(z)%*%CC + t(z)%*%DD%*%z
    if (LL){ # negative log likelihood
        detEw <- determinant(Ew,logarithm=TRUE)$modulus
        out <- (2*ns*log(2*pi) + detEw + SS)/2
    } else { # sum of squares
        out <- SS
    }
    as.numeric(out)
} # end of LL_lud2d
