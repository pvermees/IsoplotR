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
profile_weightedmean_helper <- function(sigma,X,sX,LLmax,cutoff){
    LL <- LL.sigma(sigma,X,sX)
    abs(LLmax-LL-cutoff/2)
}
# Equation 18 of Galbraith and Roberts (2012)
LL.sigma <- function(sigma,X,sX){
    wi <- 1/(sigma^2+sX^2)
    mu <- sum(wi*X)/sum(wi)
    sum(log(wi) - wi*(X-mu)^2)/2
}

ci_regression <- function(fit,model,alpha=0.05,
                          disp=TRUE,i1='b',i2='a'){
    out <- fit
    if (fit$model==1) out$fact <- nfact(alpha)
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
ci_isochron <- function(fit,model,alpha=0.05,disp=TRUE){
    ci_regression(fit,model,alpha=alpha,disp=disp,i1='age',i2='y0')
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
