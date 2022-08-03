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

maintit <- function(x,sx,n=NULL,ntit=paste0('(n=',n,')'),sigdig=2,
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
            out <- substitute(p~a%+-%b~'|'~c*u~n,lst)
        }
    } else {
        if (oerr>3){
            out <- substitute(p~a*u%+-%b*'%'~n,lst)
        } else {
            out <- substitute(p~a%+-%b*u~n,lst)
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

get.ntit <- function(x,...){ UseMethod("get.ntit",x) }
get.ntit.default <- function(x,...){
    ns <- length(x)
    nisnan <- length(which(is.na(x)))
    out <- '(n='
    if (nisnan>0) out <- paste0(out,ns-nisnan,'/')
    paste0(out,ns,')')
}
get.ntit.fissiontracks <- function(x,...){
    if (x$format<2){
        out <- get.ntit.default(x$x[,'Ns'])
    } else {
        out <- get.ntit.default(x$Ns)
    }
    out    
}
ntit.valid <- function(valid,...){
    paste0('(',sum(valid),'/',length(valid),')')
}
