regression <- function(xyz,model=1,type='york',omit=NULL,wtype='slope'){
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
        out <- list()
        fit <- stats::lm(xyz[,'Y'] ~ xyz[,'X'])
        E <- stats::vcov(fit)
        out$df <- fit$df.residual
        out$a <- c(stats::coef(fit)[1],sqrt(E[1,1]))
        out$b <- c(stats::coef(fit)[2],sqrt(E[2,2]))
        names(out$a) <- c('a','s[a]')
        names(out$b) <- c('b','s[b]')
        out$cov.ab <- E[1,2]
    } else if (identical(type,'titterington')){
        out <- list()
        fit <- stats::lm(xyz[,c('Y','Z')] ~ xyz[,'X'])
        out$df <- fit$df.residual
        out$par <- c(stats::coef(fit))
        out$cov <- stats::vcov(fit)
        parnames <- c('a','b','A','B')
        names(out$par) <- parnames
        colnames(out$cov) <- parnames
        rownames(out$cov) <- parnames
    } else {
        stop('invalid output type for model 2 regression')
    }
    out
}

model3regression <- function(xyz,type='york',wtype='slope'){
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
        lower <- init - 5
        upper <- init + 5
        fit <- stats::optim(init,LL.york,method='L-BFGS-B',
                            signs=signs,lower=lower,upper=upper,
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
        w <- stats::optimise(LL.titterington,interval=c(0,3*init),
                             xy=xyz,maximum=TRUE)$maximum
        H <- stats::optimHess(w,LL.titterington,xy=xyz)
        sw <- sqrt(solve(-H))
        out$disp <- c('w'=w,'s[w]'=sw)
        dd <- augment_titterington_errors(xyz,out$disp['w'])
        out <- c(out,titterington(dd))                
    } else {
        stop('invalid output type for model 3 regression')
    }
    out
}

york2DE <- function(XY,a,b,w=0,wtype='slope'){
    out <- list()
    ns <- nrow(XY)
    X <- XY[,'X']
    Y <- XY[,'Y']
    sX <- XY[,'sX']
    sY <- XY[,'sY']
    rXY <- XY[,'rXY']
    P <- get.york.xy(XY,a=a,b=b)
    x <- P[,'x']
    y <- P[,'y']
    E <- matrix(0,3*ns,3*ns)
    ix <- 1:ns
    iy <- (ns+1):(2*ns)
    iw <- (2*ns+1):(3*ns)
    diag(E)[ix] <- sX^2
    diag(E)[iy] <- sY^2
    diag(E)[iw] <- w^2
    E[ix,iy] <- E[iy,ix] <- diag(rXY*sX*sY)
    Jw <- matrix(0,2*ns,3*ns)
    Jw[ix,ix] <- Jw[iy,iy] <- diag(ns)
    if (wtype%in%c(0,'intercept','a')){
        Jw[ix,iw] <- -diag(P[,'dxda'])
        Jw[iy,iw] <- -diag(P[,'dyda'])
    } else {
        Jw[ix,iw] <- -diag(P[,'dxdb'])
        Jw[iy,iw] <- -diag(P[,'dydb'])
    }
    out$D <- c(X-x,Y-y)
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
LL.titterington <- function(w,xyz){
    out <- 0
    D <- augment_titterington_errors(xyz,w)
    fit <- titterington(D)
    dat <- matrix2covlist(D)
    X <- matrix(0,1,3)
    ns <- nrow(D)
    for (i in 1:ns){
        XYZ <- dat$XYZ[[i]]
        O <- dat$omega[[i]]
        E <- solve(O)
        a <- fit$par[1]
        b <- fit$par[2]
        A <- fit$par[3]
        B <- fit$par[4]
        abg <- alpha.beta.gamma(a,b,A,B,XYZ,O)
        X[1,1] <- -abg[2]/abg[1]
        X[1,2] <- XYZ[2] - a - b*XYZ[1]
        X[1,3] <- XYZ[3] - A - B*XYZ[1]
        SS <- X %*% O %*% t(X)
        detE <- determinant(E,logarithm=TRUE)$modulus
        out <- out - 1.5*log(2*pi) - 0.5*detE - 0.5*SS
    }
    out
}

augment_titterington_errors <- function(xyz,w){
    out <- xyz
    out[,'sY'] <- sqrt(xyz[,'sY']^2 + w^2)
    out[,'rXY'] <- xyz[,'rXY']*xyz[,'sY']/out[,'sY']
    out[,'rYZ'] <- xyz[,'rYZ']*xyz[,'sY']/out[,'sY']
    out
}
