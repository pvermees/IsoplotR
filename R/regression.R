regression <- function(xyz,model=1,type='york',omit=NULL){
    xyz2calc <- clear(xyz,omit)
    if (model==1)
        out <- model1regression(xyz2calc,type=type)
    else if (model==2)
        out <- model2regression(xyz2calc,type=type)
    else if (model%in%c(3,4))
        out <- model34regression(xyz2calc,type=type,model=model)
    else
        stop('invalid regression model')
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
        out <- tls(xyz[,c('X','Y')])
        out$df <- nrow(xyz)-2
        out$a <- c(out$par['a'],'s[a]'=unname(sqrt(out$cov['a','a'])))
        out$b <- c(out$par['b'],'s[b]'=unname(sqrt(out$cov['b','b'])))
        out$cov.ab <- out$cov['a','b']
    } else if (identical(type,'titterington')){
        out <- tls(xyz[,c('X','Y','Z')])
        out$df <- 2*nrow(xyz)-4
    } else {
        stop('invalid output type for model 2 regression')
    }
    out
}

# fixes signs and uses logs for numerical stability:
model34regression <- function(xyz,type='york',model=3){
    pilot <- model1regression(xyz,type=type)
    if (identical(type,'york')){
        err <- ifelse(model==3,pilot$a['s[a]'],pilot$b['s[b]'])
        fact <- max(1,sqrt(pilot$mswd))
        w <- log(fact*err)
        init <- c(pilot$a['a'],pilot$b['b'],'w'=unname(w))
        lower <- init - c(20*fact*pilot$a['s[a]'],20*fact*pilot$b['s[b]'],20)
        upper <- init + c(20*fact*pilot$a['s[a]'],20*fact*pilot$b['s[b]'],2)
        out <- stats::optim(init,LL.york,method='L-BFGS-B',
                            lower=lower,upper=upper,XY=xyz,model=model)
        x <- get.york.xy(XY=xyz,a=out$par[1],b=out$par[2],
                         w=exp(out$par[3]),model=model)[,'x']
        H <- stats::optimHess(par=c(out$par,x),fn=LL.york.ablwx,XY=xyz,model=model)
        out$cov <- solve(H)[1:3,1:3]
        out$a <- c('a'=unname(out$par['a']),'s[a]'=unname(sqrt(out$cov['a','a'])))
        out$b <- c('b'=unname(out$par['b']),'s[b]'=unname(sqrt(out$cov['b','b'])))
        out$cov.ab <- out$cov['a','b']
    } else if (identical(type,'titterington')){
        if (wtype%in%c('intercept',0,'a')) w <- sqrt(pilot$cov['a','a']*pilot$mswd)
        else if (wtype%in%c(1,'b')) w <- sqrt(pilot$cov['b','b']*pilot$mswd)
        else if (wtype%in%c(2,'A')) w <- sqrt(pilot$cov['A','A']*pilot$mswd)
        else if (wtype%in%c(3,'B')) w <- sqrt(pilot$cov['B','B']*pilot$mswd)
        else stop('illegal wtype')
        init <- c(pilot$par,'w'=unname(log(w)))
        spar <- sqrt(pilot$mswd*diag(pilot$cov))
        lower <- c(pilot$par - 5*spar,'w'=init['w']-5)
        upper <- c(pilot$par + 5*spar,'w'=init['w']+2)
        out <- stats::optim(init,LL.titterington,method='L-BFGS-B',
                            lower=lower,upper=upper,XYZ=xyz,
                            hessian=TRUE,wtype=wtype)
        out$cov <- inverthess(out$hessian)
    } else {
        stop('invalid output type for model 3 regression')
    }
    disp <- exp(out$par['w'])
    sdisp <- disp*sqrt(out$cov['w','w'])
    out$disp <- c('w'=unname(disp),'s[w]'=unname(sdisp))
    out
}

LL.york <- function(ablw,XY,model=3){
    a <- ablw[1]
    b <- ablw[2]
    w <- exp(ablw[3])
    x <- get.york.xy(XY=XY,a=a,b=b,w=w,model=model)[,'x']
    LL.york.ablwx(c(ablw,x),XY,model=model)
}
LL.york.ablwx <- function(ablwx,XY,model=3){
    ns <- nrow(XY)
    a <- ablwx[1]
    b <- ablwx[2]
    w <- exp(ablwx[3])
    x <- ablwx[4:(ns+3)]
    DE <- matrix(0,nrow=ns,ncol=5)
    colnames(DE) <- c('X-x','Y-y','vX','vY','sXY')
    DE[,'X-x'] <- XY[,'X']-x
    DE[,'Y-y'] <- XY[,'Y']-a-b*x
    DE[,'vX'] <- XY[,'sX']^2
    if (model==3){
        DE[,'vY'] <- XY[,'sY']^2 + w^2
    } else if (model==4){
        DE[,'vY'] <- XY[,'sY']^2 + (w*x)^2
    } else {
        DE[,'vY'] <- XY[,'sY']^2
    }
    DE[,'sXY'] <- XY[,'rXY']*XY[,'sX']*XY[,'sY']
    detE <- DE[,'vX']*DE[,'vY'] - DE[,'sXY']^2
    O11 <- DE[,'vY']/detE
    O22 <- DE[,'vX']/detE
    O12 <- -DE[,'sXY']/detE
    maha <- (O11*DE[,'X-x'] + O12*DE[,'Y-y'])*DE[,'X-x'] +
        (O12*DE[,'X-x'] + O22*DE[,'Y-y'])*DE[,'Y-y']
    sum(log(detE) + maha)
}
LL.york.ablwx.old <- function(ablwx,XY,model=3){
    ns <- nrow(XY)
    a <- ablwx[1]
    b <- ablwx[2]
    w <- exp(ablwx[3])
    x <- ablwx[4:(ns+3)]
    E <- matrix(0,2*ns,2*ns)
    ix <- 1:ns
    iy <- (ns+1):(2*ns)
    diag(E)[ix] <- XY[,'sX']^2
    diag(E)[iy] <- XY[,'sY']^2
    if (model==3){
        diag(E)[iy] <- diag(E)[iy] + w^2
    } else if (model==4){
        diag(E)[iy] <- diag(E)[iy] + (w*x)^2
    }
    E[ix,iy] <- E[iy,ix] <- diag(XY[,'rXY']*XY[,'sX']*XY[,'sY'])
    D <- c(XY[,'X']-x,XY[,'Y']-a-b*x)
    LL.norm(D,E)
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
LL.titterington <- function(abABw,XYZ,wtype='a'){
    DE <- titterington2DE(XYZ,a=abABw['a'],b=abABw['b'],
                          A=abABw['A'],B=abABw['B'],
                          w=exp(abABw['w']),wtype=wtype)
    LL.norm(DE$D,DE$E)
}
