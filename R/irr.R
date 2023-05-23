

LLABw <- function(ABw,XY){
    UVEOslope()
}

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
    LLui <- function(ui,A,B,lw,XYi,wtype=1){
        Ui <- log(XYi[1])
        Vi <- log(XYi[3])
        vi <- log(exp(A)+B*exp(ui))
        uvi <- c(ui,vi)
        UVi <- c(Ui,Vi)
        E11 <- (XYi[2]/XYi[1])^2
        E22 <- (XYi[4]/XYi[3])^2
        E12 <- XYi[5]*sqrt(E11*E22)
        if (wtype==1) disp <- exp(A)*exp(lw)/(exp(A)+B*exp(ui))
        else if (wtype==2) disp <- exp(ui)*exp(lw)/(exp(A)+B*exp(ui))
        else disp <- 0
        E22 <- E22 + disp^2
        detE <- E11*E22-E12^2
        Oi <- rbind(c(E22,-E12),c(-E12,E11))/detE
        SS <- stats::mahalanobis(x=UVi,center=uvi,cov=Oi,inverted=TRUE)
        log(detE) + SS
    }
    getui <- function(i,A,B,lw,XY,wtype=1){
        ll <- -100
        if (B<0) ul <- log(-exp(A)/B)-(1e-10)
        else ul <- 100
        ui <- optimise(LLui,interval=c(ll,ul),A=A,B=B,lw=lw,XYi=XY[i,],wtype=wtype)
    }
    LLABlw <- function(ABlw,XY){
        us <- sapply(X=1:nrow(XY),FUN=getui,
                     A=ABlw[1],B=ABlw[2],lw=ABlw[3],XY=XY)
        sum(unlist(us['objective',]))
    }
    yfit <- york(XY)
    init <- c(A=unname(log(yfit$a[1])),B=unname(yfit$b[1]),lw=0)
    if (wtype==1) init['lw'] <- yfit$a[2]/yfit$a[1]
    else if (wtype==2) init['lw'] <- yfit$b[2]/yfit$b[1]
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
