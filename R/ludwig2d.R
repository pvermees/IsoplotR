ludwig2d <- function(x,type=1,model=1,anchor=0,exterr=FALSE){
    if (model==1){
        out <- ludwig2d_helper(x=x,type=type,anchor=anchor,exterr=exterr)
    } else if (model==2){
        out <- ludwig2d_model2(x=x,type=type,anchor=anchor,exterr=exterr)
    } else if (model==3){
        out <- ludwig2d_helper(x=x,w=NULL,type=type,anchor=anchor,exterr=exterr)
    } else {
        stop('Invalid fit model.')
    }

    pnames <- c('t','a0','b0')
    if (model!=2) pnames <- c(pnames,'w')
    names(out$par) <- pnames
    rownames(out$cov) <- pnames
    colnames(out$cov) <- pnames

    lnames <- c('log(t)','log(a0)','log(b0)')
    if (model!=2) lnames <- c(lnames,'log(w)')
    names(out$logpar) <- lnames
    rownames(out$logcov) <- lnames
    colnames(out$logcov) <- lnames

    out$model <- model
    out$n <- length(x)
    out
}

ludwig2d_helper <- function(x,w=0,type=1,anchor=0,exterr=FALSE){

    out <- list(par=rep(NA,4),cov=matrix(NA,4,4),
                logpar=rep(NA,4),logcov=matrix(NA,4,4))    
    if (type%in%c(1,3)) i <- c(1,2)
    else if (type%in%c(2,4)) i <- c(1,3)
    else stop('Invalid isochron type.')
    if (is.null(w)){ # model 3 regression
        model1fit <- ludwig2d_helper(x=x,type=type,anchor=anchor)
        LL0 <- LL_lud2d_UPb(lta0w=model1fit$logpar[i],x=x,type=type)
        LLw <- LL_lud2d_UPb(lta0w=model1fit$logpar[i],x=x,w=0.01,type=type)
        if (LLw>LL0){ # zero dispersion
            out <- model1fit
            out$par[4] <- 0
            out$logpar[4] <- -Inf
            return(out)
        } else {
            init <- model1fit$logpar[i]
            winit <- initwlud2d(x=x,lta0=init)
            init <- c(init,winit)
            i <- c(i,4)
        }
    } else {
        model2fit <- ludwig2d_model2(x=x,type=type,anchor=anchor)
        init <- model2fit$logpar[i]
    }
    
    if (x$format%in%(4:6)){
        fit <- stats::optim(init,LL_lud2d_UPb,method='L-BFGS-B',
                            lower=init-2,upper=init+2,exterr=exterr,
                            w=w,x=x,type=type,hessian=TRUE)
        out$logpar[i] <- fit$par
        out$logcov[i,i] <- solve(fit$hessian)
        out$par[i] <- exp(out$logpar[i])
        out$cov[i,i] <- diag(out$par[i])%*%out$logcov[i,i]%*%diag(out$par[i])
        SS <- LL_lud2d_UPb(fit$par[1:2],x=x,type=type,exterr=exterr,LL=FALSE)
    } else {
        LL_UThPb <- function(lta0w,x,type=1,anchor=0){
            
        }
    }
    
    out$df <- ifelse(anchor[1]<1,length(x)-2,length(x)-1)
    out$mswd <- SS/out$df
    out$p.value <- as.numeric(1-stats::pchisq(SS,out$df))
    out
}
initwlud2d <- function(x=x,lta0,type=1){
    LL <- function(w,lta0,x,type=1){
        LL_lud2d_UPb(lta0,x=x,w=w,type=type)
    }
    stats::optim(0,LL,method='L-BFGS-B',lower=0,upper=exp(lta0[1]),
                 x=x,lta0=lta0,type=type)$par
}
LL_lud2d_UPb <- function(lta0w,x,tt=NULL,a0=NULL,w=0,
                         type=1,LL=TRUE,exterr=FALSE){
    ns <- length(x)
    np <- length(lta0w)
    if (np>2){
        tt <- exp(lta0w[1])
        a0 <- exp(lta0w[2])
        w <- exp(lta0w[3])
    } else if (np==2){
        if (!is.null(tt)){
            a0 <- exp(lta0w[1])
            w <- exp(lta0w[2])
        } else if (!is.null(a0)){
            tt <- exp(lta0w[1])
            w <- exp(lta0w[2])
        } else {
            tt <- exp(lta0w[1])
            a0 <- exp(lta0w[2])
        }
    } else if (np==1){
        if (!is.null(tt) & !is.null(a0)){
            w <- exp(lta0w)
        } else if (!is.null(tt) & !is.null(w)){
            a0 <- exp(lta0w)
        } else {
            tt <- exp(lta0w)
        }
    } else {
        if (is.null(tt) | is.null(a0)){
            stop("Missing tt and a0 values.")
        } else if (is.null(w)){
            w <- 0
        }
    }
    Y <- rep(NA,ns)
    Z <- rep(NA,ns)
    E <- matrix(0,2*ns+7,2*ns+7)
    J <- matrix(0,2*ns,2*ns+7)
    J[1:(2*ns),1:(2*ns)] <- diag(2*ns)
    D <- mclean(tt=tt,d=x$d,exterr=exterr)
    nc <- length(D$ThUi) # nc>1 if each aliquot has its own diseq correction
    j <- 1
    for (i in 1:ns){
        wd <- wetherill(x,i=i)
        if (nc>1) j <- i
        if (type==1){
            Y[i] <- wd$x['Pb206U238']
            Z[i] <- wd$x['Pb204U238']
            E[i,i] <- wd$cov['Pb206U238','Pb206U238']
            E[ns+i,ns+i] <- wd$cov['Pb204U238','Pb204U238']
            E[i,ns+i] <- wd$cov['Pb206U238','Pb204U238']
            E[ns+i,i] <- E[i,ns+i]
            J[i,2*ns+1] <- -D$dPb206U238dl38[j]  #dLdl38
            J[i,2*ns+3] <- -D$dPb206U238dl34[j]  #dLdl34
            J[i,2*ns+6] <- -D$dPb206U238dl30[j]  #dLdl30
            J[i,2*ns+7] <- -D$dPb206U238dl26[j]  #dLdl26
        } else {
            U <- iratio('U238U235')[1]
            Y[i] <- wd$x['Pb207U235']
            Z[i] <- wd$x['Pb204U238']*U
            E[i,i] <- wd$cov['Pb207U235','Pb207U235']
            E[ns+i,ns+i] <- wd$cov['Pb204U238','Pb204U238']*U^2
            E[i,ns+i] <- wd$cov['Pb207U235','Pb204U238']*U
            E[ns+i,i] <- E[i,ns+i]
            J[i,2*ns+2] <- -D$dPb207U235dl35[j]     #dKdl35
            J[i,2*ns+5] <- -D$dPb207U235dl31[j]     #dKdl31
        }
    }
    E[2*ns+1:7,2*ns+1:7] <- getEl()
    ED <- J%*%E%*%t(J)
    dY <- Y - ifelse(type==1,D$Pb206U238,D$Pb207U235)
    i1 <- 1:ns
    i2 <- ns+(1:ns)
    Jw <- matrix(0,2*ns,ns)
    diag(Jw[i1,i1]) <- -ifelse(type==1,D$dPb206U238dt,D$Pb207U235)
    Ew <- Jw%*%t(Jw)*w^2 + ED
    O <- blockinverse(Ew[i1,i1],Ew[i1,i2],
                      Ew[i1,i2],Ew[i2,i2],doall=TRUE)
    AA <- dY%*%O[i1,i1]%*%dY + dY%*%O[i1,i2]%*%Z +
        Z%*%O[i2,i1]%*%dY + Z%*%O[i2,i2]%*%Z
    BB <- a0*dY%*%O[i1,i1] + dY%*%O[i1,i2] +
        a0*Z%*%O[i2,i1] + Z%*%O[i2,i2]
    CC <- a0*O[i1,i1]%*%dY + a0*O[i1,i2]%*%Z +
        O[i2,i1]%*%dY + O[i2,i2]%*%Z
    DD <- O[i1,i1]*a0^2 + O[i1,i2]*a0 +
        O[i2,i1]*a0 + O[i2,i2]
    z <- as.vector(solve(DD+t(DD),t(BB)+CC))
    SS <- AA - BB%*%z - t(z)%*%CC + t(z)%*%DD%*%z
    if (LL){ # negative log likelihood
        detEw <- determinant(Ew,logarithm=TRUE)$modulus
        out <- (2*ns*log(2*pi) + detEw + SS)/2
    } else { # sum of squares
        out <- SS
    }
    as.numeric(out)
} # end of LL_UPb

ludwig2d_model2 <- function(x,type=1,anchor=0,exterr=FALSE){

    out <- list(par=rep(NA,3),cov=matrix(NA,3,3),
                logpar=rep(NA,3),logcov=matrix(NA,3,3))    
    
    if (type%in%c(1,3)) i <- c(1,2)
    else if (type%in%c(2,4)) i <- c(1,3)
    else stop('Invalid isochron type.')
    
    if (x$format%in%(4:6)){
        
        XY <- data2york(x,option=(3:4)[type])
        if (anchor[1]<1){
            fit <- stats::lm( XY[,'Y'] ~ XY[,'X'] )
            ab <- fit$coefficients
            covmat <- stats::vcov(fit)
        } else if (anchor[1]==1) {
            y0 <- 1/iratio('Pb206Pb204')[1]
            fit <- stats::lm( I(XY[,'Y']-y0) ~ 0 + XY[,'X'] )
            ab <- c(y0,fit$coefficients)
            covmat <- matrix(0,2,2)
            covmat[2,2] <- stats::vcov(fit)
        } else {
            x0 <- 1/mclean(anchor[2],d=x$d)$Pb206U238
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
        out$par[i] <- c(tt,Dd)
        out$cov[i,i] <- J%*%E%*%t(J)
        out$logpar[i] <- log(out$par[i])
        J <- diag(1/out$par[i])
        out$logcov[i,i] <- J%*%out$cov[i,i]%*%t(J)
        
    } else if (x$format%in%(7:8)){
        
        LL <- function(lta0,x,tt=NULL,a0=NULL,option=1,exterr=FALSE){
            ns <- length(x)
            if (length(lta0)>1){
                tt <- exp(lta0[1])
                a0 <- exp(lta0[2])
                df <- ns-2
            } else if (is.null(tt)){
                tt <- exp(lta0)
                df <- ns-1
            } else if (is.null(a0)){
                a0 <- exp(lta0)
                df <- ns-1
            } else {
                stop('You must provide initial values for both t and a0.')
            }
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
        
        model2init <- function(x,option=1,tt=NULL,a0=NULL){
            ti <- min(get.Pb206U238.age(x)[,1]) # first stab
            XY <- data2york(x,option=option,tt=ti)
            init <- c(0,0)
            if (option==6){
                init[1] <- min(get.Pb206U238.age(1/XY[,'X'])[,1])
            } else if (option==7){
                init[1] <- min(get.Pb207U235.age(1/XY[,'X'])[,1])
            } else {
                init[1] <- min(get.Pb208Th232.age(1/XY[,'X'])[,1])
            }
            init[2] <- stats::median(XY[,'Y'])
            log(init)
        }
        
        option <- (6:9)[type]
        init <- model2init(x,option=option)
        if (anchor[1]<1){
            fit <- stats::optim(init,fn=LL,x=x,option=option,
                                exterr=exterr,hessian=TRUE)
            out$logpar[i] <- fit$par
            out$logcov[i,i] <- solve(fit$hessian)
        } else if (anchor[1]==1){
            a0 <- ifelse(type%in%c(1,3),
                         1/iratio('Pb208Pb206')[1],
                         1/iratio('Pb208Pb207')[1])
            fit <- stats::optim(init[1],fn=LL,x=x,method='BFGS',a0=a0,
                                option=option,exterr=exterr,hessian=TRUE)
            out$logpar[i] <- c(fit$par,log(a0))
            out$logcov[i,i] <- matrix(0,2,2)
            out$logcov[i[1],i[1]] <- 1/fit$hessian
        } else {
            tt <- anchor[2]
            fit <- stats::optim(init[2],fn=LL,x=x,method='BFGS',tt=tt,
                                option=option,exterr=exterr,hessian=TRUE)
            out$logpar[i] <- c(log(tt),fit$par)
            out$logcov[i,i] <- matrix(0,2,2)
            out$logcov[i[2],i[2]] <- 1/fit$hessian
        }
        out$par[i] <- exp(out$logpar[i])
        J <- diag(out$par[i])
        out$cov[i,i] <- J %*% out$logcov[i,i] %*% t(J)
        
    } else {
        stop('2D ludwig regression is not available for this format')
    }
            
    out
}
