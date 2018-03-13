# based on equation 6.5 of Galbraith (2005)
profile_LL_central_disp_FT <- function(mu,sigma,y,m,alpha=0.05){
    LLmax <- LL.FT(sigma,mu,y,m)
    cutoff <- qchisq(1-alpha,1)
    if (abs(LL.FT(exp(-10),mu,y,m)-LLmax) < cutoff/2){
        sigmal <- 0
    } else {
        sigmal <- optimize(profile_central_FT_helper,
                       interval=c(0,sigma),
                       mu=mu,y=y,m=m,LLmax=LLmax,
                       cutoff=cutoff)$minimum
    }
    sigmauinit <- 2*sd(log(y)-log(m-y))
    if (abs(LL.FT(sigmauinit,mu,y,m)-LLmax) < cutoff/2){
        sigmau <- Inf
    } else {
        sigmau <- optimize(profile_central_FT_helper,
                           interval=c(sigma,sigmauinit),
                           mu=mu,y=y,m=m,LLmax=LLmax,
                           cutoff=cutoff)$minimum
    }    
    ll <- sigma - sigmal
    ul <- sigmau - sigma
    c(ll,ul)
}
profile_LL_central_disp_UThHe <- function(fit,x,alpha=0.05){
    LLmax <- LL.uvw(fit$w,fit$uvw,x=x,doSm=doSm(x))
    cutoff <- qchisq(1-alpha,1)
    doSm <- doSm(x)
    if (abs(LL.uvw(0,UVW=fit$uvw,x=x,doSm=doSm)-LLmax) < cutoff/2){
        wl <- 0
    } else {
        wl <- optimize(profile_UThHe_disp_helper,
                           interval=c(0,fit$w),x=x,UVW=fit$uvw,
                           doSm=doSm,LLmax=LLmax,cutoff=cutoff)$minimum
    }
    if (abs(LL.uvw(100,UVW=fit$uvw,x=x,doSm=doSm)-LLmax) < cutoff/2){
        wu <- Inf
    } else {
        wu <- optimize(profile_UThHe_disp_helper,
                           interval=c(fit$w,10),x=x,UVW=fit$uvw,
                           doSm=doSm,LLmax=LLmax,cutoff=cutoff)$minimum
    }
    ll <- fit$w - wl
    ul <- wu - fit$w
    c(ll,ul)
}
profile_LL_discordia_disp <- function(fit,x,alpha=0.05){
    w <- fit$w
    ta0b0 <- fit$par
    wrange <- get_lud_wrange(ta0b0,x)
    LLmax <- LL.lud.UPb.disp(w,x,ta0b0)
    cutoff <- qchisq(1-alpha,1)    
    if (abs(LL.lud.UPb.disp(wrange[1],x,ta0b0)-LLmax) < cutoff/2){
        wl <- 0
    } else {
        wl <- optimize(profile_discordia_helper,
                       interval=c(0,w),x=x,ta0b0=ta0b0,
                       LLmax=LLmax,cutoff=cutoff)$minimum
    }
    if (abs(LL.lud.UPb.disp(wrange[2],x,ta0b0)-LLmax) < cutoff/2){
        wu <- Inf
    } else {
        wu <- optimize(profile_discordia_helper,
                       interval=c(w,wrange[2]),x=x,ta0b0=ta0b0,
                       LLmax=LLmax,cutoff=cutoff)$minimum
    }
    ll <- w - wl
    ul <- wu - w
    c(ll,ul)
}

profile_LL_weightedmean_disp <- function(fit,X,sX,alpha=0.05){
    mu <- fit$mu[1]
    sigma <- fit$sigma
    LLmax <- LL.sigma(sigma,X,sX)
    cutoff <- qchisq(1-alpha,1)
    if (abs(LL.sigma(0,X,sX)-LLmax) < cutoff/2){
        sigmal <- 0
    } else {
        sigmal <- optimize(profile_weightedmean_helper,
                           interval=c(0,sigma),
                           X=X,sX=sX,LLmax=LLmax,
                           cutoff=cutoff)$minimum
    }
    if (abs(LL.sigma(10*sd(X),X,sX)-LLmax) < cutoff/2){
        sigmau <- Inf
    } else {
        sigmau <- optimize(profile_weightedmean_helper,
                           interval=c(sigma,10*sd(X)),
                           X=X,sX=sX,LLmax=LLmax,
                           cutoff=cutoff)$minimum
    }
    ll <- sigma - sigmal
    ul <- sigmau - sigma
    c(ll,ul)
}
profile_central_FT_helper <- function(sigma,mu,y,m,LLmax,cutoff){
    LL <- LL.FT(sigma,mu,y,m)
    abs(LLmax-LL-cutoff/2)
}
profile_weightedmean_helper <- function(sigma,X,sX,LLmax,cutoff){
    LL <- LL.sigma(sigma,X,sX)
    abs(LLmax-LL-cutoff/2)
}
profile_discordia_helper <- function(w,x,ta0b0,LLmax,cutoff){
    LL <- LL.lud.UPb.disp(w,x,ta0b0)
    abs(LLmax-LL-cutoff/2)
}
profile_UThHe_disp_helper <- function(w,x,UVW,doSm=FALSE,LLmax,cutoff){
    LL <- LL.uvw(w,UVW=UVW,x=x,doSm=doSm)
    abs(LLmax-LL-cutoff/2)
}

LL.FT <- function(sigma,mu,y,m){
    LL <- 0
    ns <- length(y)
    for (i in 1:ns){
        LL <- LL + log(integrate(pyumu,
                                 lower=mu-10*sigma,upper=mu+10*sigma,
                                 mu=mu,sigma=sigma,m=m[i],y=y[i])$value)
    }
    LL
}
pyumu <- function(B,mu,sigma,m,y){
    theta <- exp(B)/(1+exp(B))
    dbinom(y,m,theta)*dnorm(B,mean=mu,sd=sigma)
}

# Equation 18 of Galbraith and Roberts (2012)
LL.sigma <- function(sigma,X,sX){
    wi <- 1/(sigma^2+sX^2)
    mu <- sum(wi*X)/sum(wi)
    sum(log(wi) - wi*(X-mu)^2)/2
}

ci_regression <- function(fit,model=1,alpha=0.05,
                          disp=TRUE,i1='b',i2='a'){
    out <- fit
    if (fit$model==3) out$fact <- nfact(alpha)
    else out$fact <- tfact(alpha,fit$df)
    out[[i1]]['ci[t]'] <- out$fact*fit[[i1]]['s[t]']
    out[[i2]]['ci[y]'] <- out$fact*fit[[i2]]['s[y]']
    if (model==1 & disp){
        out[[i2]]['disp[y]'] <- sqrt(fit$mswd)*out[[i2]]['ci[y]']
    } else if (model==3) {
        out$w[c('ll','ul')] <-
            profile_LL_isochron_disp(fit,alpha=alpha)
    }
    out
    
}
ci_isochron <- function(fit,model=1,alpha=0.05,disp=TRUE){
    ci_regression(fit,model=model,alpha=alpha,disp=disp,i1='age',i2='y0')
}
profile_LL_isochron_disp <- function(fit,alpha=0.05){
    cutoff <- qchisq(1-alpha,1)
    d <- fit$d
    w <- fit$w['s']
    LLmax <- LL.isochron(w,d,type=fit$type)
    if (abs(LL.isochron(0,d,type=fit$type)-LLmax) < cutoff/2){
        wl <- 0
    } else {
        wl <- optimize(profile_isochron_helper,
                       interval=c(0,w),d=d,
                       LLmax=LLmax,cutoff=cutoff,
                       type=fit$type)$minimum
    }
    if (abs(LL.isochron(sd(d[,'Y']),d,type=fit$type)-LLmax) < cutoff/2){
        wu <- Inf
    } else {
        wu <- optimize(profile_isochron_helper,
                       interval=c(w,sd(d[,'Y'])),
                       d=d,LLmax=LLmax,cutoff=cutoff,
                       type=fit$type)$minimum
    }
    ll <- w - wl
    ul <- wu - w
    c(ll,ul)
}
profile_isochron_helper <- function(w,d,LLmax,cutoff,type='york'){
    LL <- LL.isochron(w,d,type=type)
    abs(LLmax-LL-cutoff/2)
}