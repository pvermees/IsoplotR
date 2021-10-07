ludwig2d <- function(x,type=1,model=1,anchor=0,exterr=FALSE){
    out <- ludwig2d_model2(x=x,type=type,anchor=anchor,exterr=exterr)
    out$model <- model
    out$n <- length(x)
    out
}

ludwig2d_model2 <- function(x,type=1,anchor=0,exterr=FALSE){
    out <- list(par=rep(NA,3),cov=matrix(NA,3,3),
                logpar=rep(NA,3),logcov=matrix(NA,3,3))
    if (type%in%c(1,3)) i <- c(1,2)
    else if (type%in%c(2,4)) i <- c(1,3)
    else stop('Invalid isochron type.')
    if (x$format%in%(4:6)){
        pnames <- c('t','64i','74i')
        lnames <- c('log(t)','log(64i)','log(74i)')
        option <- (3:4)[type]
        XY <- data2york(x,option=option)
        if (anchor[1]<1){
            fit <- stats::lm( XY[,'Y'] ~ XY[,'X'] )
            ab <- fit$coefficients
            covmat <- vcov(fit)
        } else if (anchor[1]==1) {
            y0 <- 1/anchor[2]
            fit <- stats::lm( I(XY[,'Y']-y0) ~ 0 + XY[,'X'] )
            ab <- c(y0,fit$coefficients)
            covmat <- matrix(0,2,2)
            covmat[2,2] <- vcov(fit)
        } else {
            x0 <- 1/mclean(anchor[2],d=x$d)$Pb206U238
            fit <- stats::lm( I(XY[,'Y']-0) ~ 0 + I(XY[,'X']-x0) )
            y0 <- -x0*fit$coefficients
            ab <- c(y0,fit$coefficients)
            covmat <- matrix(0,2,2)
            covmat[2,2] <- vcov(fit)
        }
        DP <- -ab[2]/ab[1]
        Dd <- 1/ab[1]
        if (type==1){
            tt <- get.Pb206U238.age(DP)[1]
            D <- mclean(tt,d=x$d)
            dtdPbU <- 1/D$dPb206U238dt
        } else {
            tt <- get.Pb207U235.age(DP)[1]
            D <- mclean(tt,d=x$d)
            dtdPbU <- 1/D$dPb207U235dt
        }
        J <- matrix(0,2,2)
        J[1,1] <- dtdPbU*ab[2]/ab[1]^2
        J[1,2] <- -dtdPbU/ab[1]
        J[2,1] <- -1/ab[1]^2
        out$par[i] <- c(tt,Dd)
        out$cov[i,i] <- J%*%covmat%*%t(J)
        out$logpar[i] <- log(out$par[i])
        J <- diag(1/out$par[i])
        out$logcov[i,i] <- J%*%out$cov[i,i]%*%t(J)
    } else if (x$format%in%(7:8)) {
        UThPbmisfit <- function(lta0,lt=NULL,la0=NULL,option=1,LL=TRUE,d=diseq()){
            if (length(lta0)>1){
                tt <- exp(lta0[1])
                a0 <- exp(lta0[2])
            } else if (is.null(lt)){
                tt <- exp(lta0)
                a0 <- exp(la0)
            } else if (is.null(a0)){
                tt <- exp(lt)
                a0 <- exp(lta0)
            } else {
                stop('You must provide initial values for both lt and a0.')
            }
            if (option==6){
                D <- mclean(tt,d=d)
                x0 <- 1/D$Pb206U238
                y0 <- 1/a0
            } else if (option==7){
                D <- mclean(tt,d=d)
                x0 <- 1/D$Pb207U235
                y0 <- 1/a0
            } else {
                x0 <- 1/age_to_Pb208Th232_ratio(tt)[1]
                y0 <- a0
            }
            XY <- data2york(x,option=option,tt=tt)
            yp <- y0*(1-XY[,'X']/x0)
            SS <- sum((yp-XY[,'Y'])^2)
            if (LL) out <- SS/2
            else out <- SS
            out
        }        
        option <- (6:9)[type]
        pnames <- c('t','68i','78i')
        lnames <- c('log(t)','log(68i)','log(78i)')
        if (anchor[1]<1){
            fit <- optim(c(0,0),fn=UThPbmisfit,option=option,d=x$d,hessian=TRUE)
            out$logpar[i] <- fit$par
            out$logcov[i,i] <- solve(fit$hessian)
        } else if (anchor[1]==1){
            la0 <- log(anchor[2])
            fit <- optim(0,fn=UThPbmisfit,method='BFGS',
                         la0=la0,option=option,d=x$d,hessian=TRUE)
            out$logpar[i] <- c(fit$par,la0)
            out$logcov[i,i] <- matrix(0,2,2)
            out$logcov[i[1],i[1]] <- 1/fit$hessian
        } else {
            lt <- log(anchor[2])
            fit <- optim(0,fn=UThPbmisfit,method='BFGS',
                         lt=lt,option=option,d=x$d,hessian=TRUE)
            out$logpar[i] <- c(lt,fit$par)
            out$logcov[i,i] <- matrix(0,2,2)
            out$logcov[i[2],i[2]] <- 1/fit$hessian
        }
        out$par[i] <- exp(out$logpar[i])
        J <- diag(out$par[i])
        out$cov[i,i] <- J %*% out$logcov[i,i] %*% t(J)
    } else {
        stop('2D ludwig regression is not available for this format')
    }
    names(out$par) <- pnames
    rownames(out$cov) <- pnames
    colnames(out$cov) <- pnames
    names(out$logpar) <- lnames
    rownames(out$logcov) <- lnames
    colnames(out$logcov) <- lnames
    out
}
