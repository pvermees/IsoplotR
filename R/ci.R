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
profile_weightedmean_helper <- function(sigma,X,sX,LLmax,cutoff){
    LL <- LL.sigma(sigma,X,sX)
    abs(LLmax-LL-cutoff/2)
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
geterr <- function(x,sx,oerr=5,dof=NULL){
    if (is.null(dof)) fact <- ntfact(alpha())
    else fact <- stats::qt(1-alpha()/2,df=dof)
    if (oerr==1) out <- sx
    else if (oerr==2) out <- 2*sx
    else if (oerr==3) out <- 100*sx/x
    else if (oerr==4) out <- 200*sx/x
    else if (oerr==5) out <- fact*sx
    else if (oerr==6) out <- 100*fact*sx/x
    else stop('Illegal oerr value')
}

# get significance level for error ellipses
oerr2alpha <- function(oerr=1){
    if (oerr%in%c(1,3)) out <- stats::pnorm(-1)*2
    else if (oerr%in%c(2,4)) out <- stats::pnorm(-2)*2
    else out <- alpha()
    out
}

agetit <- function(x,sx,n=NA,ntit=paste0('n=',n),sigdig=2,
                   oerr=5,units='Ma',prefix='age =',dof=NULL){
    xerr <- geterr(x,sx,oerr=oerr,dof=dof)
    rounded <- roundit(x,xerr,sigdig=sigdig,oerr=oerr,text=TRUE)
    dispersed <- (length(sx)>1)
    relerr <- (oerr %in% c(3,4,6))
    lst <- list(p=prefix,a=rounded[1],b=rounded[2],u=units,n=ntit)
    if (dispersed){
        lst$c <- rounded[3]
        if (relerr){
            out <- substitute(p~a~u%+-%b~'|'~c*'% ('*n*')',lst)
        } else {
            out <- substitute(p~a%+-%b~'|'~c~u~'('*n*')',lst)
        }
    } else {
        if (relerr){
            out <- substitute(p~a~u%+-%b*'% ('*n*')',lst)
        } else {
            out <- substitute(p~a%+-%b~u~'('*n*')',lst)
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
peaktit <- function(x,sx,p,sp,sigdig=2,oerr=5,unit='Ma',prefix=NULL){
    xerr <- geterr(x,sx,oerr=oerr)
    rounded.x <- roundit(x,xerr,sigdig=sigdig,oerr=oerr,text=TRUE)
    rounded.p <- roundit(100*p,100*sp,sigdig=sigdig,text=TRUE)
    relerr <- (oerr %in% c(3,4,6))
    lst <- list(p=prefix,a=rounded.x[1],b=rounded.x[2],c=rounded.p[1],u=unit)
    if (relerr){
        out <- substitute(p~a~u%+-%b*'% (prop='*c*'%)',lst)
    } else {
        out <- substitute(p~a%+-%b~u~'(prop='*c*'%)',lst)
    }
    out
}
