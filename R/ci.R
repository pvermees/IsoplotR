profile_LL_weightedmean_disp <- function(fit,X,sX,alpha=0.05){
    mu <- fit$mu[1]
    sigma <- fit$sigma
    LLmax <- LL.sigma(sigma,X,sX)
    cutoff <- qchisq(1-alpha,1)
    if (abs(LL.sigma(0,X,sX)-LLmax) < cutoff/2){
        sigmal <- 0
    } else {
        sigmal <- uniroot(profile_weightedmean_helper,
                          lower=0,upper=sigma,X=X,sX=sX,
                          LLmax=LLmax,cutoff=cutoff)$root
    }
    if (abs(LL.sigma(10*sd(X),X,sX)-LLmax) < cutoff/2){
        sigmau <- Inf
    } else {
        sigmau <- uniroot(profile_weightedmean_helper,
                          lower=sigma,upper=10*sd(X),X=X,sX=sX,
                          LLmax=LLmax,cutoff=cutoff)$root
    }
    cl <- sigma - sigmal
    cu <- sigmau - sigma
    c(cl,cu)
}
# Equation 18 of Galbraith and Roberts (2012)
LL.sigma <- function(sigma,X,sX){
    wi <- 1/(sigma^2+sX^2)
    mu <- sum(wi*X)/sum(wi)
    sum(log(wi) - wi*(X-mu)^2)/2
}

profile_weightedmean_helper <- function(sigma,X,sX,LLmax,cutoff){
    LL <- LL.sigma(sigma,X,sX)
    abs(LLmax-LL)-cutoff/2
}
