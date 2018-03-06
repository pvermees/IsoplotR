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

ci_isochron <- function(fit,model,alpha=0.05,disp=TRUE){
    out <- fit
    if (fit$model==1) out$fact <- nfact(1-alpha/2)
    else out$fact <- tfact(1-alpha/2,fit$df)
    out$y0['ci[y]'] <- out$fact*fit$y0['s[y]']
    out$age['ci[t]'] <- out$fact*fit$age['s[t]']
    if (model==1 & disp){
        out$y0['disp[y]'] <- sqrt(fit$mswd)*out$y0['ci[y]']
    } else if (model==3) {
        out$w[c('ll','ul')] <-
            profile_LL_isochron_disp(fit,alpha=alpha)
    }
    out
}
profile_LL_isochron_disp <- function(fit,alpha=0.05){
    cutoff <- qchisq(1-alpha,1)
    d <- fit$d
    w <- fit$w['s']
    LLmax <- LL.w(w,d)
    if (abs(LL.w(0,d)-LLmax) < cutoff/2){
        wl <- 0
    } else {
        wl <- optimize(profile_isochron_helper,
                       interval=c(0,w),d=d,
                       LLmax=LLmax,cutoff=cutoff)$minimum
    }
    if (abs(LL.w(sd(d[,'Y']),d)-LLmax) < cutoff/2){
        wu <- Inf
    } else {
        wu <- optimize(profile_isochron_helper,
                       interval=c(w,sd(d[,'Y'])),
                       d=d,LLmax=LLmax,cutoff=cutoff)$minimum
    }
    ll <- w - wl
    ul <- wu - w
    c(ll,ul)
}
profile_isochron_helper <- function(w,d,LLmax,cutoff){
    LL <- LL.w(w,d)
    abs(LLmax-LL-cutoff/2)
}
LL.w <- function(w,d){
    out <- 0
    D <- augment_york_errors(d,w)
    X <- matrix(0,1,2)
    fit <- york(D)
    P <- get.york.xy(D,fit$a[1],fit$b[1])
    for (i in 1:nrow(D)){
        E <- cor2cov2(D[i,'sX'],D[i,'sY'],D[i,'rXY'])
        X[1,1] <- D[i,'X']-P[i,1]
        X[1,2] <- D[i,'Y']-P[i,2]
        out <- out - 0.5*log(det(2*pi*E)) - 0.5 * X %*% solve(E) %*% t(X)
    }
    out
}
