dLLdui <- function(ui,A,B,lw,XYi,wtype=NA){
    Ui <- log(XYi[1])
    Vi <- log(XYi[3])
    vi <- log(exp(A)+B*exp(ui))
    w <- exp(lw)
    DD <- c(Ui-ui,Vi-vi)
    E11 <- (XYi[2]/XYi[1])^2
    E22 <- (XYi[4]/XYi[3])^2
    E12 <- XYi[5]*sqrt(E11*E22)
    EE <- rbind(c(E11,E12),c(E12,E22))
    duidui <- 1
    dvidui <- B*exp(ui)/(exp(A)+B*exp(ui))
    duidw <- 0
    d2uidwdui <- 0
    if (wtype%in%c(1,'a','intercept')) {
        dvidw <- exp(A)/(exp(A)+B*exp(ui))
        d2vidwdui <- -B*exp(A+ui)/(exp(A)+B*exp(ui))^2
    } else if (wtype%in%c(2,'b','slope')) {
        dvidw <- exp(ui)/(exp(A)+B*exp(ui))
        d2vidwdui <- -exp(A+ui)/(exp(A)+B*exp(ui))^2
    } else {
        dvidw <- 0
        d2vidwdui <- 0
    }
    Jw <- rbind(-duidw,-dvidw)
    Ew <- EE + Jw%*%(w^2)%*%t(Jw)
    Ow <- solve(Ew)
    dJwdui <- rbind(-d2uidwdui,-d2vidwdui)
    dDDdui <- rbind(-duidui,-dvidui)
    dEwdui <- dJwdui%*%(w^2)%*%t(Jw) + Jw%*%(w^2)%*%t(dJwdui)
    dLLdui <- sum(diag(Ow%*%dEwdui)) + t(dDDdui)%*%Ow%*%DD +
        t(DD)%*%Ow%*%dDDdui - t(DD)%*%(Ow%*%dEwdui%*%Ow)%*%DD
    return(dLLdui)
}

LLu <- function(u,A,B,lw,XY,wtype=NA){
    U <- log(XY[,1])
    V <- log(XY[,3])
    v <- log(exp(A)+B*exp(u))
    D1 <- U-u
    D2 <- V-v
    E11 <- (XY[,2]/XY[,1])^2
    E22 <- (XY[,4]/XY[,3])^2
    E12 <- XY[,5]*sqrt(E11*E22)
    if (wtype%in%c(1,'a','intercept')) disp <- exp(A+lw)/(exp(A)+B*exp(u))
    else if (wtype%in%c(2,'b','slope')) disp <- exp(u+lw)/(exp(A)+B*exp(u))
    else disp <- 0
    E22 <- E22 + disp^2
    detE <- E11*E22-E12^2
    O11 <- E22/detE
    O22 <- E11/detE
    O12 <- -E12/detE
    SS <- O11*D1^2 + 2*O12*D1*D2 + O22*D2^2
    LL <- (log(detE) + SS)/2
    return(LL)
}

getui <- function(i,A,B,lw,XY,wtype=NA){
    XYi <- XY[i,,drop=FALSE]
    dui <- 1e-1
    ll <- -10
    if (B<0) ul <- log(-exp(A)/B)-dui
    else ul <- 10
    mid <- min(ul-dui,log(XYi[1]))
    dLLdui_ll <- dLLdui(ll,A=A,B=B,lw=lw,XYi=XYi,wtype=wtype)
    dLLdui_mid <- dLLdui(mid,A=A,B=B,lw=lw,XYi=XYi,wtype=wtype)
    dLLdui_ul <- dLLdui(ul,A=A,B=B,lw=lw,XYi=XYi,wtype=wtype)
    # check if the gradient crosses zero
    haslower <- (dLLdui_ll*dLLdui_mid<0)
    hasupper <- (dLLdui_ul*dLLdui_mid<0)
    if (haslower){
        lui <- stats::uniroot(dLLdui,interval=c(ll,mid),A=A,B=B,lw=lw,XYi=XYi,
                              wtype=wtype,tol=.Machine$double.eps^0.5)$root
    }
    if (hasupper){
        uui <- stats::uniroot(dLLdui,interval=c(mid,ul),A=A,B=B,lw=lw,XYi=XYi,
                              wtype=wtype,tol=.Machine$double.eps^0.5)$root
    }
    if (haslower && hasupper){
        LLlui <- LLu(lui,A=A,B=B,lw=lw,XY=XYi,wtype=wtype)
        LLuui <- LLu(uui,A=A,B=B,lw=lw,XY=XYi,wtype=wtype)
        ui <- ifelse(LLlui<LLuui,lui,uui)
    } else if (haslower){
        ui <- lui
    } else if (hasupper){
        ui <- uui
    } else { # no uniroot solution => optimise the log-likelihood
        lfit <- stats::optimise(LLu,interval=c(ll,mid),
                                A=A,B=B,lw=lw,XY=XYi,wtype=wtype)
        ufit <- stats::optimise(LLu,interval=c(mid,ul),
                                A=A,B=B,lw=lw,XY=XYi,wtype=wtype)
        if (lfit$objective<ufit$objective){
            ui <- lfit$minimum
        } else {
            ui <- ufit$minimum
        }
    }
    return(ui)
}

LLABlw <- function(ABlw,XY,wtype=NA){
    np <- length(ABlw)
    A <- ABlw[1]
    B <- ABlw[2]
    if (np>2) lw <- ABlw[3]
    else lw <- -Inf
    us <- sapply(X=1:nrow(XY),FUN=getui,A=A,B=B,lw=lw,XY=XY,wtype=wtype)
    LL <- LLu(us,A=A,B=B,lw=lw,XY=XY,wtype=wtype)
    sum(LL)
}

irr <- function(XY,wtype=NA,
                control=list(maxit=10000,reltol=.Machine$double.eps),...){
    yfit <- york(XY)
    A <- unname(log(abs(yfit$a[1])))
    B <- unname(yfit$b[1])
    if (wtype%in%c(1,'a','intercept')){
        lw <- log(unname(yfit$b[2]))
        init <- c(A,B,lw)
    } else if (wtype%in%c(2,'b','slope')){
        lw <- log(unname(yfit$a[2]))
        init <- c(A,B,lw)
    } else {
        init <- c(A,B)
    }
    fit <- stats::optim(par=init,fn=LLABlw,XY=XY,hessian=TRUE,
                        wtype=wtype,control=control,...)
    np <- length(fit$par)
    if (np>2){
        w <- sqrt(exp(fit$par[3]))
        J <- diag(3)
        J[3,3] <- w/2
    } else {
        J <- diag(2)
    }
    J[1,1] <- exp(fit$par[1])
    J[2,2] <- 1
    E <- J %*% inverthess(fit$hessian) %*% t(J)
    out <- list()
    out$a <- c(exp(fit$par[1]),sqrt(E[1,1]))
    out$b <- c(fit$par[2],sqrt(E[2,2]))
    out$cov.ab <- E[1,2]
    if (np>2) out$disp <- c('w'=w,'s[w]'=sqrt(E[3,3]))
    out
}
