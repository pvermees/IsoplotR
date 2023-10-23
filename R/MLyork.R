# Maximum likelihood implementation of York regression
# This is more flexible than the original least squares method
MLyork <- function(yd,anchor=0,model=1,wtype='a'){
    tol <- 2*.Machine$double.eps
    if (anchor[1]==1 && length(anchor)>1){ # anchor intercept
        a <- anchor[2]
        init <- lm(I(yd[,'Y']-a) ~ 0 + yd[,'X'])$coefficients
        fit <- optimise(LL.MLyork.b,init,a=a,yd=yd,wtype=wtype,
                        lower=init/5,upper=init*5,tol=tol)
    } else if (anchor[1]==2 && length(anchor)>1){ # anchor slope
        b <- anchor[2]
        init <- mean(yd[,'Y']-b*yd[,'X'])
        fit <- optimise(LL.MLyork.a,init,b=b,yd=yd,wtype=wtype,
                        lower=init/5,upper=init*5,tol=tol)
    } else {
        ab <- unname(lm(yd[,'Y'] ~ yd[,'X'])$coefficients)
        init <- c(a=ab[1],b=ab[2])
        fit <- optim(init,LL.MLyork.ab,yd=yd,wtype=wtype,
                     control=list(reltol=tol))
    }
    fit
}

LL.MLyork.lw <- function(lw,ab,yd,wtype='a'){
    LL.MLyork(c(ab,lw),yd=yd,wtype=wtype)
}
LL.MLyork.a <- function(a,b,lw=-Inf,yd,wtype='a'){
    LL.MLyork(c(a,b,lw),yd=yd,wtype=wtype)
}
LL.MLyork.b <- function(b,a,lw=-Inf,yd,wtype='b'){
    LL.MLyork(c(a,b,lw),yd=yd,wtype=wtype)
}
LL.MLyork.ab <- function(ab,lw=-Inf,yd,wtype='b'){
    LL.MLyork(c(ab,lw),yd=yd,wtype=wtype)
}
LL.MLyork <- function(ablw,yd,wtype='a'){
    a <- ablw[1]
    b <- ablw[2]
    w <- exp(ablw[3])
    X <- yd[,'X']
    Y <- yd[,'Y']
    E11 <- yd[,'sX']^2
    E22 <- yd[,'sY']^2
    E12 <- yd[,'rXY']*yd[,'sX']*yd[,'sY']
    detE <- E11*E22 - E12^2
    O11 <- E22/detE
    O22 <- E11/detE
    O12 <- -E12/detE
    num <- -((O22*a-O22*Y-O12*X)*b+O12*a-O12*Y-O11*X)
    den <- O22*b^2+2*O12*b+O11
    x <- num/den
    maha <- (X-x)*(O11*(X-x)+O12*(Y-a-b*x)) +
        (Y-a-b*x)*(O12*(X-x)+O22*(Y-a-b*x))
    sum(log(detE) + maha)/2
}
