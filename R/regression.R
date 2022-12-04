regression <- function(xyz,model=1,type='york',omit=NULL,
                       wtype=ifelse(type=='york','b','a')){
    xyz2calc <- clear(xyz,omit)
    if (model==1) out <- model1regression(xyz2calc,type=type)
    else if (model==2) out <- model2regression(xyz2calc,type=type)
    else if (model==3) out <- model3regression(xyz2calc,type=type,wtype=wtype)
    else stop('invalid regression model')
    out$xyz <- xyz
    out$model <- model
    out$n <- nrow(xyz2calc)
    out$omit <- omit
    out
}

model1regression <- function(xyz,type='york'){
    if (identical(type,'york')){
        out <- york(xyz)
    } else if (identical(type,'titterington')){
        out <- titterington(xyz)
    } else {
        stop('invalid output type for model 1 regression')
    }
    out
}

model2regression <- function(xyz,type='york'){
    if (identical(type,'york')){
        pilot <- stats::lm(xyz[,'Y'] ~ xyz[,'X'])
        init <- pilot$coefficients
        names(init) <- c('a','b')
        out <- stats::optim(init,fn=LL.york.model2,XY=xyz,hessian=TRUE)
        hess <- hesscheck(out$hessian)
        out$cov <- solve(hess)
        out$df <- nrow(xyz)-2
        out$a <- c(out$par['a'],'s[a]'=unname(sqrt(out$cov['a','a'])))
        out$b <- c(out$par['b'],'s[b]'=unname(sqrt(out$cov['b','b'])))
        out$cov.ab <- out$cov['a','b']
    } else if (identical(type,'titterington')){
        pilot <- stats::lm(xyz[,c('Y','Z')] ~ xyz[,'X'])
        init <- as.vector(pilot$coefficients)
        names(init) <- c('a','b','A','B')
        out <- optim(init,fn=LL.titterington.model2,XYZ=xyz,hessian=TRUE)
        hess <- hesscheck(out$hessian)
        out$cov <- solve(hess)
        out$df <- 2*nrow(xyz)-4
    } else {
        stop('invalid output type for model 2 regression')
    }
    out
}

# fixes signs and uses logs for numerical stability:
model3regression <- function(xyz,type='york',
                             wtype=ifelse(type=='york','b','a')){
    pilot <- model1regression(xyz,type=type)
    if (identical(type,'york')){
        signs <- c('a'=sign(unname(pilot$a[1])),
                   'b'=sign(unname(pilot$b[1])),
                   'w'=1)
        a <- log(signs['a']*pilot$a[1])
        b <- log(signs['b']*pilot$b[1])
        wa <- pilot$a[2]*sqrt(pilot$mswd)
        wb <- pilot$b[2]*sqrt(pilot$mswd)
        w <- ifelse(wtype %in% c('intercept',0,'a'),log(wa),log(wb))
        init <- c('a'=unname(a),'b'=unname(b),'w'=unname(w))
        lower <- init-5
        upper <- init+5
        fit <- stats::optim(init,LL.york,method='L-BFGS-B',
                            lower=lower,upper=upper,signs=signs,
                            XY=xyz,wtype=wtype,hessian=TRUE)
        hess <- hesscheck(fit$hessian)
        fit$cov <- solve(hess)
        out <- exponentiate(fit,signs=signs)
        out$a <- c('a'=unname(out$par['a']),'s[a]'=unname(sqrt(out$cov['a','a'])))
        out$b <- c('b'=unname(out$par['b']),'s[b]'=unname(sqrt(out$cov['b','b'])))
        out$cov.ab <- out$cov['a','b']
        disp <- out$par['w']
        sdisp <- sqrt(out$cov['w','w'])
        out$disp <- c('w'=unname(disp),'s[w]'=unname(sdisp))
    } else if (identical(type,'titterington')){
        signs <- c('a'=sign(unname(pilot$par['a'])),
                   'b'=sign(unname(pilot$par['b'])),
                   'A'=sign(unname(pilot$par['A'])),
                   'B'=sign(unname(pilot$par['B'])),
                   'w'=1)
        a <- log(signs['a']*pilot$par['a'])
        b <- log(signs['b']*pilot$par['b'])
        A <- log(signs['A']*pilot$par['A'])
        B <- log(signs['B']*pilot$par['B'])
        if (wtype%in%c('intercept',0,'a')) w <- sqrt(pilot$cov['a','a']*pilot$mswd)
        else if (wtype%in%c(1,'b')) w <- sqrt(pilot$cov['b','b']*pilot$mswd)
        else if (wtype%in%c(2,'A')) w <- sqrt(pilot$cov['A','A']*pilot$mswd)
        else if (wtype%in%c(3,'B')) w <- sqrt(pilot$cov['B','B']*pilot$mswd)
        else stop('illegal wtype')
        init <- c('a'=unname(a),'b'=unname(b),
                  'A'=unname(A),'B'=unname(B),'w'=unname(w))
        lower <- init-5
        upper <- init+5
        fit <- stats::optim(init,LL.titterington,method='L-BFGS-B',
                            XYZ=xyz,lower=lower,upper=upper,
                            hessian=TRUE,wtype=wtype)
        hess <- hesscheck(fit$hessian)
        fit$cov <- solve(hess)
        out <- exponentiate(fit,signs=signs)
        disp <- out$par['w']
        sdisp <- sqrt(out$cov['w','w'])
        out$disp <- c('w'=unname(disp),'s[w]'=unname(sdisp))
    } else {
        stop('invalid output type for model 3 regression')
    }
    out
}

york2DE <- function(XY,a,b,w=0,wtype='slope'){
    out <- list()
    ns <- nrow(XY)
    P <- get.york.xy(XY,a=a,b=b)
    E <- matrix(0,3*ns,3*ns)
    ix <- 1:ns
    iy <- (ns+1):(2*ns)
    iw <- (2*ns+1):(3*ns)
    diag(E)[ix] <- XY[,'sX']^2
    diag(E)[iy] <- XY[,'sY']^2
    diag(E)[iw] <- w^2
    E[ix,iy] <- E[iy,ix] <- diag(XY[,'rXY']*XY[,'sX']*XY[,'sY'])
    Jw <- matrix(0,2*ns,3*ns)
    Jw[ix,ix] <- Jw[iy,iy] <- diag(ns)
    if (wtype%in%c(0,'intercept','a')){
        Jw[ix,iw] <- -diag(P[,'dxda'])
        Jw[iy,iw] <- -diag(P[,'dyda'])
    } else {
        Jw[ix,iw] <- -diag(P[,'dxdb'])
        Jw[iy,iw] <- -diag(P[,'dydb'])
    }
    out$D <- c(XY[,'X']-P[,'x'],XY[,'Y']-P[,'y'])
    out$E <- Jw %*% E %*% t(Jw)
    out
}
LL.york <- function(labw,XY,signs=labw/labw,wtype='slope'){
    a <- signs['a']*exp(labw['a'])
    b <- signs['b']*exp(labw['b'])
    w <- exp(labw['w'])
    DE <- york2DE(XY,a=a,b=b,w=w,wtype=wtype)
    LL.norm(DE$D,DE$E)
}
LL.york.model2 <- function(ab,XY){
    dem <- deming(a=ab['a'],b=ab['b'],x=XY[,'X'],y=XY[,'Y'])
    SS2LL(SS=sum(dem$d^2),nn=nrow(XY))
}

titterington2DE <- function(XYZ,a,b,A,B,w=0,wtype='a'){
    out <- list()
    ns <- nrow(XYZ)
    P <- get.titterington.xyz(XYZ,a=a,b=b,A=A,B=B)
    D <- c(XYZ[,'X']-P[,'x'],XYZ[,'Y']-P[,'y'],XYZ[,'Z']-P[,'z'])
    E <- matrix(0,4*ns,4*ns)
    ix <- 1:ns
    iy <- (ns+1):(2*ns)
    iz <- (2*ns+1):(3*ns)
    iw <- (3*ns+1):(4*ns)
    diag(E)[ix] <- XYZ[,'sX']^2
    diag(E)[iy] <- XYZ[,'sY']^2
    diag(E)[iz] <- XYZ[,'sZ']^2
    diag(E)[iw] <- w^2
    E[ix,iy] <- E[iy,ix] <- diag(XYZ[,'rXY']*XYZ[,'sX']*XYZ[,'sY'])
    E[ix,iz] <- E[iz,ix] <- diag(XYZ[,'rXZ']*XYZ[,'sX']*XYZ[,'sZ'])
    E[iy,iz] <- E[iz,iy] <- diag(XYZ[,'rYZ']*XYZ[,'sY']*XYZ[,'sZ'])
    Jw <- matrix(0,3*ns,4*ns)
    Jw[ix,ix] <- Jw[iy,iy] <- Jw[iz,iz] <- diag(ns)
    if (wtype%in%c(0,'intercept','a')){
        Jw[ix,iw] <- -diag(P[,'dxda'])
        Jw[iy,iw] <- -diag(P[,'dyda'])
        Jw[iz,iw] <- -diag(P[,'dzda'])
    } else if (wtype%in%c(1,'b')){
        Jw[ix,iw] <- -diag(P[,'dxdb'])
        Jw[iy,iw] <- -diag(P[,'dydb'])
        Jw[iz,iw] <- -diag(P[,'dzdb'])
    } else if (wtype%in%c(2,'A')){
        Jw[ix,iw] <- -diag(P[,'dxdA'])
        Jw[iy,iw] <- -diag(P[,'dydA'])
        Jw[iz,iw] <- -diag(P[,'dzdA'])
    } else if (wtype%in%c(3,'B')){
        Jw[ix,iw] <- -diag(P[,'dxdB'])
        Jw[iy,iw] <- -diag(P[,'dydB'])
        Jw[iz,iw] <- -diag(P[,'dzdB'])
    } else {
        stop('invalid wtype')
    }
    out$D <- c(XYZ[,'X']-P[,'x'],XYZ[,'Y']-P[,'y'],XYZ[,'Z']-P[,'z'])
    out$E <- Jw %*% E %*% t(Jw)
    out
}
LL.titterington <- function(labABw,XYZ,signs=labABw/labABw,wtype='a'){
    a <- signs['a']*exp(labABw['a'])
    b <- signs['b']*exp(labABw['b'])
    A <- signs['A']*exp(labABw['A'])
    B <- signs['B']*exp(labABw['B'])
    w <- exp(labABw['w'])
    DE <- titterington2DE(XYZ,a=a,b=b,A=A,B=B,w=w,wtype=wtype)
    LL.norm(DE$D,DE$E)
}
LL.titterington.model2 <- function(abAB,XYZ){
    demXY <- deming(a=abAB['a'],b=abAB['b'],x=XYZ[,'X'],y=XYZ[,'Y'])
    demXZ <- deming(a=abAB['A'],b=abAB['B'],x=XYZ[,'X'],y=XYZ[,'Z'])
    dd <- cbind(demXY$d,demXZ$d)
    E <- stats::cov(dd)
    sum(apply(dd,1,LL.norm,covmat=E))
}
