# Maximum likelihood implementation of York regression
# This is more flexible than the original least squares method
MLyork <- function(yd,anchor=0,model=1,wtype='a',
                   tol=.Machine$double.eps^0.5,
                   control=list(reltol=tol)){
    p <- c(NA,NA,-Inf)
    E <- matrix(0,3,3)
    names(p) <- rownames(E) <- colnames(E) <- c('a','b','lw')
    ns <- nrow(yd)
    if (anchor[1]==1 & length(anchor)>1){ # anchor intercept
        p['a'] <- anchor[2]
        lmfit <- stats::lm(I(yd[,'Y']-p['a']) ~ 0 + yd[,'X'])
        if (model==2){
            fit <- tls(yd[,c('X','Y')],anchor=anchor)
            p['b'] <- fit$par[2]
            E['b','b'] <- fit$cov[2,2]
        } else if (model==3 & (length(anchor)<3 || anchor[3]<=0)){
            binit <- unname(lmfit$coefficients)
            lwinit <- log(stats::vcov(lmfit))/2
            init <- c(b=binit,lw=lwinit)
            lower <- c(min(binit*c(2,1/2)),lwinit-2)
            upper <- c(max(binit*c(2,1/2)),lwinit+2)
            fit <- contingencyfit(init,LL.MLyork.blw,a=p['a'],
                                  yd=yd,lower=lower,upper=upper)
            pnames <- c('b','lw')
            p[pnames] <- fit$par
            E[pnames,pnames] <- inverthess(fit$hessian)
            df <- ns-2
        } else {
            binit <- unname(lmfit$coefficients)
            interval <- sort(binit*c(2,1/2))
            lw <- ifelse(model==3 & length(anchor)>2, log(anchor[3]), -Inf)
            fit <- stats::optimise(LL.MLyork.b,a=p['a'],lw=lw,
                                   yd=yd,interval=interval,tol=tol)
            p['b'] <- fit$minimum
            p['lw'] <- lw
            H <- stats::optimHess(fit$minimum,LL.MLyork.b,
                                  a=p['a'],lw=lw,yd=yd,control=control)
            E['b','b'] <- inverthess(H)
            E['lw','lw'] <- 0
            df <- ns-1
        }
    } else if (anchor[1]==2 & length(anchor)>1){ # anchor slope
        p['b'] <- anchor[2]
        init <- unname(stats::median(yd[,'Y']-p['b']*yd[,'X']))
        if (model==2){
            fit <- tls(yd[,c('X','Y')],anchor=anchor)
            p['a'] <- fit$par[1]
            E['a','a'] <- fit$cov[1,1]
        } else if (model==3 & (length(anchor)<3 || anchor[3]<=0)){
            ainit <- unname(lmfit$coefficients)
            lwinit <- log(stats::vcov(lmfit))/2
            init <- c(a=ainit,lw=lwinit)
            lower <- c(min(ainit*c(2,1/2)),lwinit-2)
            upper <- c(max(ainit*c(2,1/2)),lwinit+2)
            fit <- contingencyfit(init,LL.MLyork.alw,b=p['b'],
                                  yd=yd,lower=lower,upper=upper)
            pnames <- c('a','lw')
            p[pnames] <- fit$par
            E[pnames,pnames] <- inverthess(fit$hessian)
            df <- ns-2
        } else {
            interval <- sort(c(init/2,init*2))
            lw <- ifelse(model==3 & length(anchor)>2, log(anchor[3]), -Inf)
            fit <- stats::optimise(LL.MLyork.a,init,b=p['b'],lw=lw,
                                   yd=yd,interval=interval,tol=tol)
            p['a'] <- fit$minimum
            p['lw'] <- lw
            H <- stats::optimHess(fit$minimum,LL.MLyork.a,
                                  b=p['b'],lw=lw,yd=yd,control=control)
            E['a','a'] <- inverthess(H)
            E['lw','lw'] <- 0
            df <- ns-1
        }
    } else {
        lmfit <- stats::lm(yd[,'Y'] ~ yd[,'X'])
        ab <- unname(lmfit$coefficients)
        if (model==2){
            i <- c('a','b')
            fit <- tls(yd[,c('X','Y')])
            p[i] <- fit$par
            E[i,i] <- fit$cov
        } else if (model==3){
            init <- c(a=ab[1],b=ab[2],lw=0)
            fit <- stats::optim(init,LL.MLyork.ablw,yd=yd,wtype=wtype,
                                hessian=TRUE,control=control)
            p <- fit$par
            E <- inverthess(fit$hessian)
        } else { # this is equivalent to ordinary York regression
            init <- c(a=ab[1],b=ab[2])
            fit <- stats::optim(init,LL.MLyork.ab,yd=yd,
                                hessian=TRUE,control=control)
            p[1:2] <- fit$par
            E <- inverthess(fit$hessian)
            SS <- LL.MLyork.ab(fit$par,yd=yd,returnval='SS')
            df <- ns-2
        }
    }
    if (model==1){
        SS <- LL.MLyork.ablw(p,yd=yd,returnval='SS')
        fit <- append(fit,getMSWD(X2=SS,df=df))
    }
    fit$model <- model
    fit$n <- ns
    fit$par <- p
    fit$cov <- E
    fit$a <- c(p['a'],'s[a]'=unname(sqrt(E['a','a'])))
    fit$b <- c(p['b'],'s[b]'=unname(sqrt(E['b','b'])))
    fit$cov.ab <- E['a','b']
    if (model==3){
        fit$w <- c('w'=unname(exp(p['lw'])),
                   's[w]'=unname(exp(p['lw'])*sqrt(E['lw','lw'])))
    }
    fit
}

LL.MLyork.a <- function(a,b,lw=-Inf,yd,returnval='LL'){
    LL.MLyork.ablw(c(a,b,lw),yd=yd,returnval=returnval)
}
LL.MLyork.b <- function(b,a,lw=-Inf,yd,returnval='LL'){
    LL.MLyork.ablw(c(a,b,lw),yd=yd,returnval=returnval)
}
LL.MLyork.ab <- function(ab,lw=-Inf,yd,returnval='LL'){
    LL.MLyork.ablw(c(ab,lw),yd=yd,returnval=returnval)
}
LL.MLyork.alw <- function(alw,b,yd,returnval='LL'){
    LL.MLyork.ablw(c(alw[1],b,alw[2]),yd=yd,
                   wtype='b',returnval=returnval)
}
LL.MLyork.blw <- function(blw,a,yd,returnval='LL'){
    LL.MLyork.ablw(c(a,blw),yd=yd,
                   wtype='a',returnval=returnval)
}
LL.MLyork.ablw <- function(ablw,yd,wtype='a',returnval='LL'){
    a <- ablw[1]
    b <- ablw[2]
    w <- exp(ablw[3])
    if (wtype%in%c(2,'slope','b')){
        misfit <- function(x,yd,a,b,w,SS=FALSE){
            E <- MLY.getE(yd=yd,w=w,wtype=wtype,x=x)
            O <- MLY.getO(E)
            S <- MLY.maha(a=a,b=b,x=x,X=yd[,'X'],Y=yd[,'Y'],O=O)
            if (SS) return(sum(S))
            detE <- E[,1]*E[,2]-E[,3]^2
            sum(log(detE)+S)/2
        }
        fit <- stats::optim(yd[,'X'],misfit,yd=yd,a=a,b=b,w=w)
        if (returnval=='LL'){
            out <- fit$value
        } else {
            x <- fit$par
            if (returnval=='x') out <- x
            else if (returnval=='SS') out <- misfit(x=x,yd=yd,a=a,b=b,w=w,SS=TRUE)
            else stop("Invalid returnval.")
        }
    } else { # wtype%in%c(1,'intercept','a')
        X <- yd[,'X']
        Y <- yd[,'Y']
        E <- MLY.getE(yd=yd,w=w,wtype=wtype)
        O <- MLY.getO(E)
        O11 <- O[,1]; O22 <- O[,2]; O12 <- O[,3]
        num <- -((O22*a-O22*Y-O12*X)*b+O12*a-O12*Y-O11*X)
        den <- O22*b^2+2*O12*b+O11
        x <- num/den
        if (returnval=='x'){
            out <- x
        } else {
            Ew <- MLY.getE(yd,w=w,wtype=wtype,x=x)
            S <- MLY.maha(a=a,b=b,x=x,X=yd[,'X'],Y=yd[,'Y'],O=MLY.getO(Ew))
            if (returnval=='SS'){
                out <- sum(S)
            } else {
                detEw <- Ew[,1]*Ew[,2]-Ew[,3]^2
                out <- sum(log(detEw)+S)/2
            }
        }
    }
}

MLY.getE <- function(yd,w=0,wtype='a',x=0){
    E11 <- yd[,'sX']^2
    E22 <- yd[,'sY']^2
    E12 <- yd[,'rXY']*yd[,'sX']*yd[,'sY']
    if (wtype%in%c(2,'slope','b')) E22 <- E22 + (w*x)^2
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
