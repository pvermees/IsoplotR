LRisochron <- function(x,...){ UseMethod("LRisochron",x) }
LRisochron.default <- function(x,left=TRUE,...){
    lower <- c(-Inf,0,0,min(x[,'Y']))
    upper <- c(Inf,1,Inf,max(x[,'Y']))
    y0i <- x[which.min(x[,'X']),'Y']
    gi <- 0
    propi <- 0.5
    sigi <- sd((y0i-x[,'Y'])/x[,'X'])
    init <- c(gi,propi,sigi,y0i)
    fun <- ifelse(left,yd2ratios.left,yd2ratios.right)
    fit <- contingencyfit(par=init,fn=get_LRisochron_L,
                          lower=lower,upper=upper,
                          yd=x,fun=fun,hessian=TRUE)
    covmat <- inverthess(fit$hessian)
    a <- fit$par[4]
    sa <- sqrt(covmat[4,4])
    b <- ifelse(left,-1,1) * fit$par[1]
    sb <- sqrt(covmat[1,1])
    list(a=c('a'=a,'s[a]'=sa),
         b=c('b'=b,'s[b]'=sb),
         cov.ab=covmat[1,4],
         model=1)
}
LRisochron.PbPb <- function(x,inverse=TRUE,anchor=0,...){
    yd <- data2york(x,inverse=TRUE)
    gi <- 0
    propi <- 0.5
    y0i <- yd[which.min(yd[,'X']),'Y']
    sigi <- sd((y0i-yd[,'Y'])/yd[,'X'])
    x1 <- y1 <- y0 <- NULL
    if (anchor[1]==1 & length(anchor)>1){
        y0 <- age2ratio(tt=anchor[2],ratio="Pb207Pb206")[1]
        init <- c(gi,propi,sigi)
        lower <- c(-Inf,0,0)
        upper <- c(Inf,1,Inf)
    } else if (anchor[1]==2){
        Pb74 <- iratio('Pb207Pb204')[1]
        Pb64 <- iratio('Pb206Pb204')[1]
        y1 <- Pb74/Pb64
        x1 <- 1/Pb64
        init <- c(propi,sigi,y0i)
        lower <- c(0,0,min(yd[,'Y']))
        upper <- c(1,Inf,max(yd[,'Y']))
    } else {
        init <- c(gi,propi,sigi,y0i)
        lower <- c(-Inf,0,0,min(yd[,'Y']))
        upper <- c(Inf,1,Inf,max(yd[,'Y']))
    }
    fit <- contingencyfit(par=init,fn=get_LRisochron_L,
                          lower=lower,upper=upper,
                          fun=yd2ratios.left,
                          yd=yd,y0=y0,x1=x1,y1=y1,
                          hessian=TRUE)
    covmat <- inverthess(fit$hessian)
    np <- length(fit$par)
    if (is.null(y0)){
        a <- fit$par[np]
        sa <- sqrt(covmat[np,np])
    } else {
        a <- y0
        sa <- 0
        cov.ab <- 0
    }
    if (is.null(x1) | is.null(y1)){
        b <- -fit$par[1]
        sb <- sqrt(covmat[1,1])
        if (is.null(y0)) cov.ab <- covmat[1,np]
    } else {
        b <- (y1-a)/x1
        sb <- sa/x1
        rho <- -1
        cov.ab <- rho*sa*sb
    }
    out <- list(a=c('a'=unname(a),'s[a]'=unname(sa)),
                b=c('b'=unname(b),'s[b]'=unname(sb)),
                cov.ab=unname(cov.ab))
    if (inverse){
        return(out)
    } else {
        return(invertfit(out,type="d"))
    }
}

LRisochron.ThU <- function(x,type=1,UTh=NULL,...){
    if (x$format<3){
        stop("Rightmost isochrons are only available for ThU formats 3 and 4.")
    }
    yd <- data2york(x,type=1)
    lower <- c(-Inf,0,0)
    upper <- c(Inf,1,Inf)
    gi <- 0
    propi <- 0.5
    y0i <- ifelse(is.null(UTh),1,1/UTh)
    sigi <- sd((yd[,'Y']-y0i)/(yd[,'X']-y0i))
    init <- c(gi,propi,sigi)
    if (is.null(UTh)){
        lower <- append(lower,min(yd[,'Y']))
        upper <- append(lower,max(yd[,'Y']))
        init <- append(init,y0i)
        y0 <- NULL
    } else {
        y0 <- y0i
    }
    fit <- contingencyfit(par=init,fn=get_LRisochron_L,
                          lower=lower,upper=upper,
                          yd=yd,y0=y0,fun=yd2ratios.ThU,
                          hessian=TRUE)
    covmat <- inverthess(fit$hessian)
    b <- fit$par[1]
    sb <- sqrt(covmat[1,1])
    if (is.null(UTh)){
        y0 <- fit$par[4]
        sy0 <- sqrt(covmat[4,4])
        E <- errorprop(J11=-y0,J12=1-b,
                       J21=1,J22=0,
                       E11=covmat[1,1],
                       E22=covmat[4,4],
                       E12=covmat[1,4])
    } else {
        sy0 <- 0
        E <- errorprop(J11=-y0,J12=1-b,
                       J21=1,J22=0,
                       E11=covmat[1,1],
                       E22=0,E12=0)
    }
    a <- y0*(1-b)
    sa <- sqrt(E[1,1])
    cov.ab <- E[1,2]
    out <- list(y0=c('y0'=unname(y0),'s[y0]'=unname(sy0)),
                a=c('a'=unname(a),'s[a]'=unname(sa)),
                b=c('b'=unname(b),'s[b]'=unname(sb)),
                cov.ab=unname(cov.ab))
    if (type==1){
        return(out)
    } else {
        return(invertfit(out,type="p"))
    }
}

yd2ratios.left <- function(yd,y0){
    r <- (y0-yd[,'Y'])/yd[,'X']
    sr <- sqrt(errorprop1x2(J1=-r/yd[,'X'],
                            J2=-1/yd[,'X'],
                            E11=yd[,'sX']^2,
                            E22=yd[,'sY']^2,
                            E12=yd[,'rXY']*yd[,'sX']*yd[,'sY']))
    cbind(r,sr)
}
yd2ratios.right <- function(yd,y0){
    r <- (yd[,'Y']-y0)/yd[,'X']
    sr <- sqrt(errorprop1x2(J1=-r/yd[,'X'],
                            J2=1/yd[,'X'],
                            E11=yd[,'sX']^2,
                            E22=yd[,'sY']^2,
                            E12=yd[,'rXY']*yd[,'sX']*yd[,'sY']))
    cbind(r,sr)    
}
yd2ratios.ThU <- function(yd,y0){
    r <- (yd[,'Y']-y0)/(yd[,'X']-y0)
    sr <- sqrt(errorprop1x2(J1=-r/(yd[,'X']-y0),
                            J2=-y0/(yd[,'X']-y0),
                            E11=yd[,'sX']^2,
                            E22=yd[,'sY']^2,
                            E12=yd[,'rXY']*yd[,'sX']*yd[,'sY']))
    cbind(r,sr)
}

get_LRisochron_L <- function(pars,yd,y0=NULL,x1=NULL,y1=NULL,fun=yd2ratios){
    np <- length(pars)
    if (is.null(y0)){
        y0 <- pars[np]
    }
    if (is.null(x1) || is.null(y1)){
        gam <- pars[1]
    } else {
        gam <- (y0-y1)/x1
    }
    prop <- pars[2]
    sig <- pars[3]
    mu <- gam
    zs <- fun(yd=yd,y0=y0)
    z <- zs[,1]
    s <- zs[,2]
    AA  <- prop/sqrt(2*pi*s^2)
    BB <- -0.5*((z-gam)/s)^2
    CC <- (1-prop)/sqrt(2*pi*(sig^2+s^2))
    mu0 <- (mu/sig^2 + z/s^2)/(1/sig^2 + 1/s^2)
    s0 <- 1/sqrt(1/sig^2 + 1/s^2)
    DD <- 1-stats::pnorm((gam-mu0)/s0)
    EE <- 1-stats::pnorm((gam-mu)/sig)
    FF <- -0.5*((z-mu)^2)/(sig^2+s^2)
    fu <- AA*exp(BB) + CC*(DD/EE)*exp(FF)
    fu[fu<.Machine$double.xmin] <- .Machine$double.xmin
    fu[fu>.Machine$double.xmax] <- .Machine$double.xmax
    sum(-log(fu))
}

