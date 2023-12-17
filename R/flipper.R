# helper function for isochron() that flips X- and Y- axis and back
# if necessary to facilitate anchored regression and model-3 fits
flipper <- function(x,inverse=FALSE,hide=NULL,omit=NULL,model=1,
                    wtype=0,anchor=0,type='p',y0rat,t2DPfun,...){
    flip <- FALSE
    xyz <- data2york(x,inverse=inverse)
    if (model<3 & anchor[1]<1){
        fitinverse <- inverse
        d2calc <- clear(xyz,hide,omit)
        fit <- regression(d2calc,model=model)
    } else if (anchor[1]==1){
        wtype <- 1 # override
        fitinverse <- FALSE
        d2calc <- invertandclean(x=x,inverse=inverse,
                                 fitinverse=fitinverse,
                                 hide=hide,omit=omit)
        if (!missing(y0rat)) anchor[2:3] <- iratio(y0rat)
        if (model>1) fit <- MLyork(d2calc,anchor=anchor,model=model)
        else fit <- anchoredYork(d2calc,y0=anchor[2],sy0=anchor[3])
    } else if (anchor[1]==2){
        wtype <- 2 # override
        fitinverse <- TRUE
        d2calc <- invertandclean(x=x,inverse=inverse,fitinverse=fitinverse,
                                 hide=hide,omit=omit)
        flip <- (type=='p')
        if (flip) d2calc[,c('X','sX','Y','sY','rXY')] <- d2calc[,c(3,4,1,2,5)]
        if (missing(t2DPfun)){
            DP <- anchor[2:3]
        } else {
            st <- ifelse(length(anchor)<3,0,anchor[3])
            DP <- do.call(t2DPfun,args=list(t=anchor[2],st=st,...))
        }
        if (type=='d') y0 <- DP
        else y0 <- quotient(X=DP[1],sX=DP[2],Y=1,sY=0,sXY=0)
        if (model>1){
            anchor <- c(1,y0)
            fit <- MLyork(d2calc,anchor=anchor,model=model)
        } else {
            fit <- anchoredYork(d2calc,y0=y0[1],sy0=y0[2])
        }
    } else if (wtype==1){
        fitinverse <- FALSE
        d2calc <- invertandclean(x=x,inverse=inverse,
                                 fitinverse=fitinverse,
                                 hide=hide,omit=omit)
        fit <- MLyork(d2calc,model=model,wtype='a')
    } else if (wtype==2){
        fitinverse <- TRUE
        d2calc <- invertandclean(x=x,inverse=inverse,
                                 fitinverse=fitinverse,
                                 hide=hide,omit=omit)
        flip <- (type=='p')
        if (flip) yd[,c('X','sX','Y','sY','rXY')] <- yd[,c(3,4,1,2,5)]
        fit <- MLyork(d2calc,model=model,wtype='a')
    } else {
        stop("Invalid anchor and/or wtype value.")
    }
    fit$anchor <- anchor
    out <- list()
    if (flip){
        out$flippedfit <- fit
        fit <- unflipfit(fit)
    }
    inverted <- (inverse != fitinverse)
    if (inverted){
        out$invertedfit <- fit
        fit <- invertfit(fit,type=type,wtype=wtype)
    }
    out <- append(out,fit)
    out$xyz <- xyz
    out
}

invertandclean <- function(x,inverse,fitinverse,hide,omit){
    if (is.other(x) & inverse!=fitinverse){
        yd <- normal2inverse(data2york(x))
    } else {
        yd <- data2york(x,inverse=fitinverse)
    }
    clear(yd,hide,omit)
}

# the purpose of flipping and inverting is always to use the intercept
# to improve the fit, therefore the inverse operations always attribute
# any overdispersion to the incoming intercept
unflipfit <- function(fit){
    out <- fit
    a <- -fit$a[1]/fit$b[1]
    b <- 1/fit$b[1]
    J11 <- -1/fit$b[1]
    J12 <- -a/fit$b[1]
    J21 <- 0
    J22 <- -b/fit$b[1]
    E11 <- fit$a[2]^2
    E22 <- fit$b[2]^2
    E12 <- fit$cov.ab
    vcovab <- errorprop(J11,J12,J21,J22,E11,E22,E12)
    out$a <- c(a,sqrt(vcovab[1]))
    out$b <- c(b,sqrt(vcovab[2]))
    out$cov.ab <- unname(vcovab[3])
    out
}
invertfit <- function(fit,type="p",wtype=0){
    out <- fit
    if (type%in%c(1,"p")){
        a <- 1/fit$a[1]
        b <- -fit$b[1]/fit$a[1]
        J11 <- -a/fit$a[1]
        J12 <- 0
        J21 <- -b/fit$a[1]
        J22 <- -1/fit$a[1]
        E11 <- fit$a[2]^2
        E22 <- fit$b[2]^2
        E12 <- fit$cov.ab
        vcovab <- errorprop(J11,J12,J21,J22,E11,E22,E12)
        out$a <- c(a,sqrt(vcovab[1]))
        out$b <- c(b,sqrt(vcovab[2]))
        out$cov.ab <- vcovab[3]
    } else if (type%in%c(2,"d")){
        out$a <- fit$b
        out$b <- fit$a
    } else {
        stop("Invalid isochron type.")
    }
    names(out$a) <- c('a','s[a]')
    names(out$b) <- c('b','s[b]')
    out
}

anchoredYork <- function(x,y0=0,sy0=0){
    eps <- .Machine$double.eps
    X <- rbind(x,c(0,eps,y0,max(sy0,eps),0))
    out <- yorkhelper(X,np=1)
    if (y0==0) out$a[1] <- 0
    if (sy0==0) out$a[2] <- 0
    out$model <- 1
    out$n <- nrow(x)
    out
}
