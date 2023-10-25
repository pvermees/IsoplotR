# Maximum likelihood implementation of York regression
# This is more flexible than the original least squares method
MLyork <- function(yd,anchor=0,model=1,wtype='a'){
    tol <- 2*.Machine$double.eps
    p <- c(NA,NA,-Inf)
    E <- matrix(0,3,3)
    names(p) <- rownames(E) <- colnames(E) <- c('a','b','lw')
    ns <- nrow(yd)
    if (anchor[1]==1 && length(anchor)>1){ # anchor intercept
        p['a'] <- anchor[2]
        init <- lm(I(yd[,'Y']-p['a']) ~ 0 + yd[,'X'])$coefficients
        if (model==2){
            i <- 'b'
            fit <- tls(yd[,c('X','Y')],anchor=anchor)
            p[i] <- fit$par[2]
            E[i,i] <- fit$cov[2,2]
        } else if (model==3){ # wtype = 2
            i <- c('b','lw')
            init <- c(b=init,lw=0)
            fit <- optim(init,LL.MLyork.blw,a=p['a'],yd=yd,wtype=wtype,
                         control=list(reltol=tol),hessian=TRUE)
            p[i] <- fit$par
            E[i,i] <- inverthess(fit$hessian)
        } else {
            i <- 'b'
            fit <- optimise(LL.MLyork.b,a=p['a'],yd=yd,
                            lower=init/5,upper=init*5,tol=tol)
            p[i] <- fit$minimum
            H <- optimHess(fit$minimum,LL.MLyork.b,a=p['a'],yd=yd)
            E[i,i] <- inverthess(H)
            df <- ns-1
        }
    } else if (anchor[1]==2 && length(anchor)>1){ # anchor slope
        p['b'] <- anchor[2]
        init <- mean(yd[,'Y']-p['b']*yd[,'X'])
        if (model==2){
            i <- 'a'
            fit <- tls(yd[,c('X','Y')],anchor=anchor)
            p[i] <- fit$par[1]
            E[i,i] <- fit$cov[1,1]
        } else if (model==3){ # wtype = 1
            i <- c('a','lw')
            init <- c(a=init,lw=0)
            fit <- optim(init,LL.MLyork.alw,b=p['b'],yd=yd,wtype=wtype,
                         control=list(reltol=tol),hessian=TRUE)
            p[i] <- fit$par
            E[i,i] <- inverthess(fit$hessian)
        } else {
            i <- 'a'
            fit <- optimise(LL.MLyork.a,b=p['b'],yd=yd,
                            lower=init/5,upper=init*5,tol=tol)
            p[i] <- fit$minimum
            H <- optimHess(fit$minimum,LL.MLyork.a,b=p['b'],yd=yd)
            E[i,i] <- inverthess(H)
            df <- ns-1
        }
    } else {
        ab <- unname(lm(yd[,'Y'] ~ yd[,'X'])$coefficients)
        if (model==2){
            i <- c('a','b')
            fit <- tls(yd[,c('X','Y')])
            p[i] <- fit$par
            E[i,i] <- fit$cov
        } else if (model==3){
            init <- c(a=ab[1],b=ab[2],lw=0)
            fit <- optim(init,LL.MLyork.ablw,yd=yd,wtype=wtype,
                         control=list(reltol=tol),hessian=TRUE)
            p <- fit$par
            E <- inverthess(fit$hessian)
        } else { # this is equivalent to ordinary York regression
            init <- c(a=ab[1],b=ab[2])
            fit <- optim(init,LL.MLyork.ab,yd=yd,
                         control=list(reltol=tol),hessian=TRUE)
            p <- fit$par
            E <- inverthess(fit$hessian)
            SS <- LL.MLyork.ab(fit$par,yd=yd,SS=TRUE)
            df <- ns-2
        }
    }
    if (model==1){
        SS <- LL.MLyork.ablw(p,yd=yd,SS=TRUE)
        fit <- append(fit,getMSWD(X2=SS,df=df))
    }
    fit$par <- p
    fit$cov <- E
    fit$a <- c(p['a'],'s[a]'=unname(sqrt(E['a','a'])))
    fit$b <- c(p['b'],'s[b]'=unname(sqrt(E['b','b'])))
    fit$cov.ab <- E['a','b']
    fit$disp <- c('w'=unname(exp(p['lw'])),
                  's[w]'=unname(exp(p['lw'])*sqrt(E['lw','lw'])))
    fit
}

LL.MLyork.a <- function(a,b,yd){
    LL.MLyork.ablw(c(a,b,-Inf),yd=yd)
}
LL.MLyork.b <- function(b,a,yd){
    LL.MLyork.ablw(c(a,b,-Inf),yd=yd)
}
LL.MLyork.ab <- function(ab,yd){
    LL.MLyork.ablw(c(ab,-Inf),yd=yd)
}
LL.MLyork.alw <- function(alw,b,yd,wtype='a'){
    if (wtype%in%c(2,'slope','b')){
        stop("It doesn't make sense to attribute ",
             "overdispersion to a fixed slope.")
    }
    LL.MLyork.ablw(c(alw[1],b,alw[2]),yd=yd,wtype=wtype)
}
LL.MLyork.blw <- function(blw,a,yd,wtype='b'){
    if (wtype%in%c(1,'intercept','a')){
        stop("It doesn't make sense to attribute ",
             "overdispersion to a fixed intercept.")
    }
    LL.MLyork.ablw(c(a,blw),yd=yd,wtype=wtype)
}

LL.MLyork.ablw <- function(ablw,yd,wtype='a',SS=FALSE){
    a <- ablw[1]
    b <- ablw[2]
    w <- exp(ablw[3])
    x <- MLY.getx(yd,a=a,b=b,w=w,wtype=wtype)
    E <- MLY.getE(yd,w=w,wtype=wtype,x=x)
    detE <- E[,1]*E[,2]-E[,3]^2
    S <- MLY.maha(a=a,b=b,x=x,X=yd[,'X'],Y=yd[,'Y'],O=MLY.getO(E))
    if (SS) return(sum(S))
    else return(sum(log(detE)+S)/2)
}

MLY.getE <- function(yd,w=0,wtype='a',x=NA){
    E11 <- yd[,'sX']^2
    E22 <- yd[,'sY']^2
    E12 <- yd[,'rXY']*yd[,'sX']*yd[,'sY']
    if (wtype%in%c('intercept',1,'b')) E22 <- E22 + (w*x)^2
    else E22 <- E22 + w^2
    cbind(E11,E22,E12)
}
MLY.getO <- function(E){
    invertcovmat(vx=E[,1],vy=E[,2],sxy=E[,3])
}
MLY.maha <- function(yd,a,b,x,X,Y,O){
    O11 <- O[,1]; O22 <- O[,2]; O12 <- O[,3]
    (X-x)*(O11*(X-x)+O12*(Y-a-b*x)) + (Y-a-b*x)*(O12*(X-x)+O22*(Y-a-b*x))
}
MLY.getx <- function(yd,a,b,w=0,wtype='a'){
    if (wtype%in%c('intercept',1,'b')){
        misfit <- function(x,yd,a,b,w){
            E <- MLY.getE(yd=yd,w=w,wtype='b',x=x)
            O <- MLY.getO(E)
            SS <- MLY.maha(a=a,b=b,x=x,X=yd[,'X'],Y=yd[,'Y'],O=O)
            detE <- E[,1]*E[,2]-E[,3]^2
            sum(log(detE)+SS)/2
        }
        x <- optim(yd[,'X'],misfit,yd=yd,a=a,b=b,w=w)$par
    } else {
        X <- yd[,'X']
        Y <- yd[,'Y']
        E <- MLY.getE(yd=yd,w=w,wtype=wtype)
        O <- MLY.getO(E)
        O11 <- O[,1]; O22 <- O[,2]; O12 <- O[,3]
        num <- -((O22*a-O22*Y-O12*X)*b+O12*a-O12*Y-O11*X)
        den <- O22*b^2+2*O12*b+O11
        x <- num/den
    }
    x
}
