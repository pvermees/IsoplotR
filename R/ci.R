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
    sigmauinit <- 2*stats::sd(log(y+0.5)-log(m-y))
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
profile_UThHe_disp_helper <- function(w,x,UVW,doSm=FALSE,LLmax,cutoff){
    LL <- LL.uvw(w,UVW=UVW,x=x,doSm=doSm)
    abs(LLmax-LL-cutoff/2)
}

LL.FT <- function(par,y,m){
    mu <- par[1]
    sigma <- exp(par[2])
    LL <- 0
    ns <- length(y)
    for (i in 1:ns){
        LL <- LL + log(stats::integrate(pyumu,lower=mu-5*sigma,
                                        upper=mu+5*sigma,mu=mu,
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

ci_regression <- function(fit,i1='b',i2='a'){
    out <- fit
    out[[i1]]['ci[t]'] <- ntfact(fit$alpha)*fit[[i1]]['s[t]']
    out[[i2]]['ci[y]'] <- ntfact(fit$alpha)*fit[[i2]]['s[y]']
    if (inflate(fit)){
        out[[i2]]['disp[y]'] <- ntfact(fit$alpha,fit)*fit[[i2]]['s[y]']
    } else if (fit$model==3) {
        out$w[c('ll','ul')] <- profile_LL_isochron_disp(fit)
    }
    out    
}
ci_isochron <- function(fit){
    ci_regression(fit,i1='age',i2='y0')
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
    LL <- LL.isochron(stats::sd(xyz[,'Y']),xyz=xyz,type=fit$type)
    if (abs(LL-LLmax) < cutoff/2){
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

ci_log2lin_lud <- function(fit,fact=1){
    lx <- fit$logpar['log(w)']
    slx <- sqrt(fit$logcov['log(w)','log(w)'])
    if (is.finite(lx)){
        ll <- exp(lx - fact*slx)
        ul <- exp(lx + fact*slx)
    } else {
        ll <- 0
        ul <- NA
    }
    c(exp(lx),ll,ul)
}

ntfact <- function(alpha=0.05,mswd=NULL,df=NULL){
    if (is.null(mswd)){
        if (is.null(df)){
            out <- stats::qnorm(1-alpha/2)
        } else if (df>0){
            out <- stats::qt(1-alpha/2,df=df)
        }
    } else if (mswd$df>0){
        out <- stats::qt(1-alpha/2,df=mswd$df)*sqrt(mswd$mswd)
    } else {
        out <- stats::qnorm(1-alpha/2)
    }
    out
}

inflate <- function(fit){
    if (is.null(fit$model)) fit$model <- 1
    (fit$model==1) && (fit$p.value<alpha())
}

# df != NULL for fission track data
geterr <- function(x,sx,oerr=5,df=NULL){
    if (is.null(df)) fact <- ntfact(alpha())
    else fact <- stats::qt(1-alpha()/2,df=df)
    if (oerr==1) out <- sx
    else if (oerr==2) out <- 2*sx
    else if (oerr==3) out <- 100*sx/x
    else if (oerr==4) out <- 200*sx/x
    else if (oerr==5) out <- fact*sx
    else if (oerr==6) out <- 100*fact*sx/x
    else stop('Illegal oerr value')
}

oerr2alpha <- function(oerr=1){
    if (oerr%in%c(1,3)) out <- stats::pnorm(-1)*2
    else if (oerr%in%c(2,4)) out <- stats::pnorm(-2)*2
    else out <- alpha()
    out
}

agetit <- function(x,sx,n=NA,ntit=paste0('n=',n),sigdig=2,
                   oerr=5,prefix='age =',df=NULL){
    xerr <- geterr(x,sx,oerr=oerr,df=df)
    rounded <- roundit(x,xerr,sigdig=sigdig,oerr=oerr,text=TRUE)
    dispersed <- (length(sx)>1)
    relerr <- (oerr %in% c(3,4,6))
    lst <- list(p=prefix,a=rounded[1],b=rounded[2],u=units,n=ntit)
    if (dispersed){
        lst$c <- rounded[3]
        if (relerr){
            out <- substitute(p~a~'Ma'%+-%b~'|'~c*'% ('*n*')',lst)
        } else {
            out <- substitute(p~a%+-%b~'|'~c~'Ma ('*n*')',lst)
        }
    } else {
        if (relerr){
            out <- substitute(p~a~'Ma'%+-%b*'% ('*n*')',lst)
        } else {
            out <- substitute(p~a%+-%b~'Ma ('*n*')',lst)
        }
    }
    out
}
mswdtit <- function(mswd,p,sigdig=2){
    substitute('MSWD ='~a~', p('*chi^2*') ='~b,
               list(a=signif(mswd,sigdig),b=signif(p,sigdig)))
}
disptit <- function(w,sw,sigdig=2,oerr=5,prefix='dispersion ='){
    werr <- geterr(w,sw,oerr=oerr)
    rounded <- roundit(w,werr,sigdig=sigdig,oerr=oerr,text=TRUE)
    lst <- list(p=prefix,a=rounded[1],b=rounded[2])
    relerr <- (oerr %in% c(3,4,6))
    if (relerr){
        out <- substitute(p~a%+-%b*'%',lst)
    } else if (werr/w<0.5){
        out <- substitute(p~a%+-%b,lst)
    } else {
        lst$b <- signif(w-exp(log(w)-werr/w),sigdig)
        lst$c <- signif(exp(log(w)+werr/w)-w,sigdig)
        out <- substitute(p~a+b-c,lst)
    }
    out
}
