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
geterr <- function(x,sx,oerr=3,dof=NULL,absolute=FALSE){
    if (oerr>3 & absolute) oerr <- oerr-3
    if (is.null(dof)) fact <- ntfact(alpha())
    else fact <- stats::qt(1-alpha()/2,df=dof)
    if (oerr==1) out <- sx
    else if (oerr==2) out <- 2*sx
    else if (oerr==3) out <- fact*sx
    else if (oerr==4) out <- abs(100*sx/x)
    else if (oerr==5) out <- abs(200*sx/x)
    else if (oerr==6) out <- abs(100*fact*sx/x)
    else stop('Illegal oerr value')
    out
}

# get significance level for error ellipses
oerr2alpha <- function(oerr=1){
    if (oerr%in%c(1,4)) out <- stats::pnorm(-1)*2
    else if (oerr%in%c(2,5)) out <- stats::pnorm(-2)*2
    else out <- alpha()
    out
}

agetit <- function(x,sx,n=NULL,ntit=paste0('(n=',n,')'),sigdig=2,
                   oerr=3,units=' Ma',prefix='age =',dof=NULL){
    xerr <- geterr(x,sx,oerr=oerr,dof=dof)
    rounded <- roundit(x,xerr,sigdig=sigdig,oerr=oerr,text=TRUE)
    dispersed <- (length(sx)>1)
    lst <- list(p=prefix,a=rounded[1],b=rounded[2],u=units,n=ntit)
    if (dispersed){
        lst$c <- rounded[3]
        if (oerr>3){
            out <- substitute(p~a*u%+-%b~'|'~c*'%'~n,lst)
        } else {
            out <- substitute(p~a%+-%b~'|'~c~u~n,lst)
        }
    } else {
        if (oerr>3){
            out <- substitute(p~a*u%+-%b*'%'~n,lst)
        } else {
            out <- substitute(p~a%+-%b~u~n,lst)
        }
    }
    out
}
mswdtit <- function(mswd,p,sigdig=2){
    substitute('MSWD ='~a*', p('*chi^2*') ='~b,
               list(a=signif(mswd,sigdig),b=signif(p,sigdig)))
}
disptit <- function(w,sw,sigdig=2,oerr=3,units='',prefix='dispersion ='){
    if (w>0){
        werr <- geterr(w,sw,oerr=oerr)
        rounded <- roundit(w,werr,sigdig=sigdig,oerr=oerr,text=TRUE)
    } else {
        w <- 0
        werr <- NA
        rounded <- c(w,werr)
    }
    lst <- list(p=prefix,a=rounded[1],b=rounded[2],u=units)
    if (oerr>3){
        out <- substitute(p~a~u%+-%b*'%',lst)
    } else if (is.na(werr) | sw/w<0.5){
        out <- substitute(p~a%+-%b~u,lst)
    } else {
        lst$b <- signif(exp(log(w)+werr/w)-w,sigdig)
        lst$c <- signif(w-exp(log(w)-werr/w),sigdig)
        out <- substitute(p~a+b-c~u,lst)
    }
    out
}
peaktit <- function(x,sx,p,sigdig=2,oerr=3,unit='Ma',prefix=NULL){
    xerr <- geterr(x,sx,oerr=oerr)
    rounded.x <- roundit(x,xerr,sigdig=sigdig,oerr=oerr,text=TRUE)
    rounded.p <- signif(100*p,sigdig)
    lst <- list(p=prefix,a=rounded.x[1],b=rounded.x[2],c=rounded.p[1],u=unit)
    if (oerr>3){
        out <- substitute(p~a~u%+-%b*'% (prop='*c*'%)',lst)
    } else {
        out <- substitute(p~a%+-%b~u~'(prop='*c*'%)',lst)
    }
    out
}
