# based on equation 6.5 of Galbraith (2005)
profile_LL_central_disp_FT <- function(mu,sigma,y,m,alpha=0.05){
    LLmax <- LL.FT(mu=mu,sigma=sigma,y=y,m=m)
    cutoff <- stats::qchisq(1-alpha,1)
    if (abs(LL.FT(mu=mu,sigma=exp(-10),y=y,m=m)-LLmax) < cutoff/2){
        sigmal <- 0
    } else {
        sigmal <- stats::optimize(profile_central_FT_helper,
                                  interval=c(0,sigma),
                                  mu=mu,y=y,m=m,LLmax=LLmax,
                                  cutoff=cutoff)$minimum
    }
    sigmauinit <- 2*stats::sd(log(y)-log(m-y))
    if (abs(LL.FT(mu=mu,sigma=sigmauinit,y=y,m=m)-LLmax) < cutoff/2){
        sigmau <- Inf
    } else {
        sigmau <- stats::optimize(profile_central_FT_helper,
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
    cutoff <- stats::qchisq(1-alpha,1)
    doSm <- doSm(x)
    if (abs(LL.uvw(0,UVW=fit$uvw,x=x,doSm=doSm)-LLmax) < cutoff/2){
        wl <- 0
    } else {
        wl <- stats::optimize(profile_UThHe_disp_helper,
                              interval=c(0,fit$w),x=x,UVW=fit$uvw,
                              doSm=doSm,LLmax=LLmax,cutoff=cutoff)$minimum
    }
    if (abs(LL.uvw(100,UVW=fit$uvw,x=x,doSm=doSm)-LLmax) < cutoff/2){
        wu <- Inf
    } else {
        wu <- stats::optimize(profile_UThHe_disp_helper,
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
    wrange <- c(0,1)
    LLmax <- LL.lud.disp(w=w,x=x,ta0b0=ta0b0)
    cutoff <- stats::qchisq(1-alpha,1)    
    if (abs(LL.lud.disp(w=wrange[1],x=x,ta0b0=ta0b0)-LLmax) < cutoff/2){
        wl <- 0
    } else {
        wl <- stats::optimize(profile_discordia_helper,interval=c(wrange[1],w),
                              x=x,ta0b0=ta0b0,LLmax=LLmax,cutoff=cutoff)$minimum
    }
    if (abs(LL.lud.disp(w=wrange[2],x=x,ta0b0=ta0b0)-LLmax) < cutoff/2){
        wu <- Inf
    } else {
        wu <- stats::optimize(profile_discordia_helper,interval=c(w,wrange[2]),
                              x=x,ta0b0=ta0b0,LLmax=LLmax,cutoff=cutoff)$minimum
    }
    ll <- w - wl
    ul <- wu - w
    c(ll,ul)
}

profile_LL_weightedmean_disp <- function(fit,X,sX,alpha=0.05){
    mu <- fit$mu[1]
    sigma <- fit$sigma
    LLmax <- LL.sigma(sigma,X,sX)
    cutoff <- stats::qchisq(1-alpha,1)
    if (abs(LL.sigma(0,X,sX)-LLmax) < cutoff/2){
        sigmal <- 0
    } else {
        sigmal <- stats::optimize(profile_weightedmean_helper,
                                  interval=c(0,sigma),
                                  X=X,sX=sX,LLmax=LLmax,
                                  cutoff=cutoff)$minimum
    }
    if (abs(LL.sigma(10*stats::sd(X),X,sX)-LLmax) < cutoff/2){
        sigmau <- Inf
    } else {
        sigmau <- stats::optimize(profile_weightedmean_helper,
                                  interval=c(sigma,10*stats::sd(X)),
                                  X=X,sX=sX,LLmax=LLmax,
                                  cutoff=cutoff)$minimum
    }
    ll <- sigma - sigmal
    ul <- sigmau - sigma
    c(ll,ul)
}
profile_central_FT_helper <- function(sigma,mu,y,m,LLmax,cutoff){
    # update mu from sigma using section 3.9.1 from Galbraith (2005)
    theta <- 1/(1+exp(-mu))
    wj <- m/(theta*(1-theta) + (m-1)*(theta^2)*((1-theta)^2)*(sigma^2))
    theta <- sum(wj*y/m)/sum(wj)
    mu <- log(theta)-log(1-theta)
    LL <- LL.FT(mu=mu,sigma=sigma,y=y,m=m)
    abs(LLmax-LL-cutoff/2)
}
profile_weightedmean_helper <- function(sigma,X,sX,LLmax,cutoff){
    LL <- LL.sigma(sigma,X,sX)
    abs(LLmax-LL-cutoff/2)
}
profile_discordia_helper <- function(w,x,ta0b0,LLmax,cutoff){
    LL <- LL.lud.disp(w=w,x=x,ta0b0=ta0b0)
    abs(LLmax-LL-cutoff/2)
}
profile_UThHe_disp_helper <- function(w,x,UVW,doSm=FALSE,LLmax,cutoff){
    LL <- LL.uvw(w,UVW=UVW,x=x,doSm=doSm)
    abs(LLmax-LL-cutoff/2)
}

LL.FT <- function(mu,sigma,y,m){
    LL <- 0
    ns <- length(y)
    for (i in 1:ns){
        LL <- LL + log(stats::integrate(pyumu,lower=mu-10*sigma,
                                        upper=mu+10*sigma,mu=mu,
                                        sigma=sigma,m=m[i],y=y[i])$value)
    }
    LL
}
pyumu <- function(B,mu,sigma,m,y){
    theta <- exp(B)/(1+exp(B))
    stats::dbinom(y,m,theta)*stats::dnorm(B,mean=mu,sd=sigma)
}

# Equation 18 of Galbraith and Roberts (2012)
LL.sigma <- function(sigma,X,sX){
    wi <- 1/(sigma^2+sX^2)
    mu <- sum(wi*X)/sum(wi)
    sum(log(wi) - wi*(X-mu)^2)/2
}

ci_regression <- function(fit,disp=TRUE,i1='b',i2='a'){
    out <- fit
    if (fit$model==3) out$fact <- nfact(fit$alpha)
    else out$fact <- tfact(fit$alpha,fit$df)
    out[[i1]]['ci[t]'] <- out$fact*fit[[i1]]['s[t]']
    out[[i2]]['ci[y]'] <- out$fact*fit[[i2]]['s[y]']
    if (fit$model==1 & disp){
        out[[i2]]['disp[y]'] <- sqrt(fit$mswd)*out[[i2]]['ci[y]']
    } else if (fit$model==3) {
        out$w[c('ll','ul')] <- profile_LL_isochron_disp(fit)
    }
    out
    
}
ci_isochron <- function(fit,disp=TRUE){
    ci_regression(fit,disp=disp,i1='age',i2='y0')
}
profile_LL_isochron_disp <- function(fit){
    cutoff <- stats::qchisq(1-fit$alpha,1)
    w <- fit$w['s']
    xyz <- fit$xyz
    LLmax <- LL.isochron(w,xyz,type=fit$type)
    if (abs(LL.isochron(0,xyz,type=fit$type)-LLmax) < cutoff/2){
        wl <- 0
    } else {
        wl <- stats::optimize(profile_isochron_helper,interval=c(0,w),xyz=xyz,
                              LLmax=LLmax,cutoff=cutoff,type=fit$type)$minimum
    }
    if (abs(LL.isochron(stats::sd(xyz[,'Y']),xyz=xyz,type=fit$type)-LLmax) < cutoff/2){
        wu <- Inf
    } else {
        wu <- stats::optimize(profile_isochron_helper,
                              interval=c(w,stats::sd(xyz[,'Y'])),
                              xyz=xyz,LLmax=LLmax,cutoff=cutoff,
                              type=fit$type)$minimum
    }
    ll <- w - wl
    ul <- wu - w
    c(ll,ul)
}
profile_isochron_helper <- function(w,xyz,LLmax,cutoff,type='york'){
    LL <- LL.isochron(w,xyz,type=type)
    abs(LLmax-LL-cutoff/2)
}
