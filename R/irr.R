irr.model1 <- function(XY){
    getUVEOslope <- function(XY,A=0,B=0){
        ns <- nrow(XY)
        U <- log(XY[,1])
        V <- log(XY[,3])
        E11 <- (XY[,2]/XY[,1])^2
        E22 <- (XY[,4]/XY[,3])^2
        E12 <- XY[,5]*sqrt(E11*E22)
        detE <- E11*E22-E12^2
        O11 <- E22/detE
        O22 <- E11/detE
        O12 <- -E12/detE
        lam <- (E11+E22)/2 + sqrt(4*E12^2+(E11-E22)^2)/2
        slope <- (lam-E11)/E12
        out <- cbind(U,V,E11,E22,E12,detE,O11,O22,O12,slope)
        colnames(out) <- c('U','V','E11','E22','E12','detE',
                           'O11','O22','O12','slope')
        out
    }
    limits <- function(i=1,A,B,UVEOslope){
        ns <- nrow(UVEOslope)
        U <- UVEOslope[i,'U']
        V <- UVEOslope[i,'V']
        b <- UVEOslope[i,'slope']
        a <- V - b*U
        out <- -100
        if ((B*b)>0){
            x1 <- (B*exp(-A)/b)^(1/(b-1))
            y1 <- x1^b
            y2 <- exp(A-a) + B*exp(-a)*x1
        }
        if (B<0){
            ul <- log(-exp(A)/B)-(1e-10)
            x2 <- -exp(-a)/B
            if (b<0 && x2<x1 && y2>y1) out <- append(out,U)
        } else {
            ul <- 100
            if (b>0 && b<1 && y1>y2) out <- append(out,U)
        }
        out <- append(out,ul)
        out
    }
    SSui <- function(ui,A,B,UVi,Oi){
        vi <- log(exp(A)+B*exp(ui))
        uvi <- c(ui,vi)
        stats::mahalanobis(x=UVi,center=uvi,cov=Oi,inverted=TRUE)
    }
    # first derivative of SSui (quicker if there is just one intersection)
    uimisfit <- function(ui,A=A,B=B,UVi,Oi){
        evi <- exp(A)+B*exp(ui)
        vi <- log(evi)
        Di <- UVi-c(ui,vi)
        dDidui <- -c(1,B*exp(ui)/evi)
        Di %*% Oi %*% dDidui + dDidui %*% Oi %*% Di
    }
    getui <- function(i,A,B,UVEOslope){
        lims <- limits(i=i,A=A,B=B,UVEOslope=UVEOslope)
        UVi <- UVEOslope[i,c('U','V')]
        Oi <- rbind(c(UVEOslope[i,'O11'],UVEOslope[i,'O12']),
                    c(UVEOslope[i,'O12'],UVEOslope[i,'O22']))
        if (length(lims)<3){
            fit <- uniroot(uimisfit,A=A,B=B,UVi=UVi,Oi=Oi,
                           lower=lims[1],upper=lims[2])
            SS <- SSui(fit$root,A,B,UVi,Oi)
            u <- fit$root
        } else {
            lfit <- optimise(SSui,A=A,B=B,UVi=UVi,Oi=Oi,
                             lower=lims[1],upper=lims[2])
            ufit <- optimise(SSui,A=A,B=B,UVi=UVi,Oi=Oi,
                             lower=lims[2],upper=lims[3])
            if (lfit$objective<ufit$objective){
                SS <- lfit$objective
                u <- lfit$minimum
            } else {
                SS <- ufit$objective
                u <- ufit$minimum
            }
        }
        out <- c(u=u,SS=SS)
    }
    SSAB <- function(AB,UVEOslope){
        us <- sapply(X=1:nrow(UVEOslope),FUN=getui,
                     A=AB[1],B=AB[2],UVEOslope=UVEOslope)
        sum(us['SS',])
    }
    yfit <- york(XY)
    UVEOslope <- getUVEOslope(XY)
    init <- c(log(yfit$a[1]),yfit$b[1])
    optim(par=init,fn=SSAB,UVEOslope=UVEOslope,hessian=TRUE)
}

irr.model3 <- function(XY,wtype=1){
    dLLdui <- function(ui,A,B,lw,XYi,wtype=1){
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
        if (wtype==1) {
            dvidw <- exp(A)/(exp(A)+B*exp(ui))
            d2vidwdui <- -B*exp(A+ui)/(exp(A)+B*exp(ui))^2
        } else if (wtype==2) {
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
    LLu <- function(u,A,B,lw,XY,wtype=1){
        U <- log(XY[,1])
        V <- log(XY[,3])
        v <- log(exp(A)+B*exp(u))
        D1 <- U-u
        D2 <- V-v
        E11 <- (XY[,2]/XY[,1])^2
        E22 <- (XY[,4]/XY[,3])^2
        E12 <- XY[,5]*sqrt(E11*E22)
        if (wtype==1) disp <- exp(A+lw)/(exp(A)+B*exp(u))
        else if (wtype==2) disp <- exp(u+lw)/(exp(A)+B*exp(u))
        else disp <- 0
        E22 <- E22 + disp^2
        detE <- E11*E22-E12^2
        O11 <- E22/detE
        O22 <- E11/detE
        O12 <- -E12/detE
        SS <- O11*D1^2 + 2*O12*D1*D2 + O22*D2^2
        (log(detE) + SS)/2
    }
    getui <- function(i,A,B,lw,XY,wtype=1){
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
            lui <- uniroot(dLLdui,interval=c(ll,mid),A=A,B=B,lw=lw,
                           XYi=XYi,wtype=wtype)$root
        }
        if (hasupper){
            uui <- uniroot(dLLdui,interval=c(mid,ul),A=A,B=B,lw=lw,
                           XYi=XYi,wtype=wtype)$root
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
            lfit <- optimise(LLu,interval=c(ll,mid),
                             A=A,B=B,lw=lw,XY=XYi,wtype=wtype)
            ufit <- optimise(LLu,interval=c(mid,ul),
                             A=A,B=B,lw=lw,XY=XYi,wtype=wtype)
            if (lfit$objective<ufit$objective){
                ui <- lfit$minimum
            } else {
                ui <- ufit$minimum
            }
        }
        return(ui)
    }
    LLABlw <- function(ABlw,XY){
        A <- ABlw[1]
        B <- ABlw[2]
        lw <- ABlw[3]
        us <- sapply(X=1:nrow(XY),FUN=getui,A=A,B=B,lw=lw,XY=XY)
        LL <- LLu(us,A=A,B=B,lw=lw,XY=XY,wtype=wtype)
        sum(LL)
    }
    yfit <- york(XY)
    init <- c(A=unname(log(yfit$a[1])),B=unname(yfit$b[1]),lw=0)
    if (wtype==1) init['lw'] <- log(yfit$a[2])
    else if (wtype==2) init['lw'] <- log(yfit$b[2])
    else init['lw'] <- -Inf    
    fit <- optim(par=init,fn=LLABlw,XY=XY,hessian=TRUE)
    fit
}

irr <- function(XY,wtype=0){
    if (wtype==0){
        fit <- irr.model1(XY)
    } else {
        fit <- irr.model3(XY,wtype=wtype)
    }
    J <- 0*fit$hessian
    J[1,1] <- exp(fit$par[1])
    J[2,2] <- 1
    if (wtype>0) {
        w <- sqrt(exp(fit$par[3]))
        J[3,3] <- w/2
    }
    E <- J %*% inverthess(fit$hessian) %*% t(J)
    out <- list()
    out$a <- c(exp(fit$par[1]),sqrt(E[1,1]))
    out$b <- c(fit$par[2],sqrt(E[2,2]))
    out$cov.ab <- E[1,2]
    if (wtype>0) out$w <- c(w,sqrt(E[3,3]))
    out
}
