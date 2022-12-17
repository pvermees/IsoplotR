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

#' @title Confidence intervals
#' @description Given a parameter estimate and its standard error,
#'     calculate the corresponding 1-sigma, 2-sigma or
#'     \eqn{100(1-\alpha)\%} confidence interval, in absolute or
#'     relative units.
#' @param x scalar estimate
#' @param sx scalar or vector with the standard error of x without and
#'     (optionally) with \eqn{\sqrt{MSWD}} overdispersion multiplier.
#' @param oerr indicates if the confidence intervals should be
#'     reported as:
#' 
#' \code{1}: 1\eqn{\sigma} absolute uncertainties.
#' 
#' \code{2}: 2\eqn{\sigma} absolute uncertainties.
#' 
#' \code{3}: absolute (1-\eqn{\alpha})\% confidence intervals, where
#' \eqn{\alpha} equales the value that is stored in
#' \code{settings('alpha')}.
#'
#' \code{4}: 1\eqn{\sigma} relative uncertainties (\eqn{\%}).
#' 
#' \code{5}: 2\eqn{\sigma} relative uncertainties (\eqn{\%}).
#'
#' \code{6}: relative (1-\eqn{\alpha})\% confidence intervals, where
#' \eqn{\alpha} equales the value that is stored in
#' \code{settings('alpha')}.
#' 
#' @param df (optional) number of degrees of freedom. Only used if
#'     \code{sx} is a vector.
#' @param absolute logical. Returns absolute uncertainties even if
#'     \code{oerr} is greater than 3. Used for some internal
#'     \code{IsoplotR} functions.
#' @return A scalar or vector of the same size as \code{sx}.
#' @details Several of \code{IsoplotR}'s plotting functions (including
#'     \code{\link{isochron}}, \code{\link{weightedmean}},
#'     \code{\link{concordia}}, \code{\link{radialplot}} and
#'     \code{\link{helioplot}}) return lists of parameters and their
#'     standard errors. For `model-1' fits, if the data pass a
#'     Chi-square test of homogeneity, then just one estimate for the
#'     standard error is reported.  This estimate can be converted to
#'     a confidence interval by multiplication with the appropriate
#'     quantile of a Normal distribution. Datasets that fail the
#'     Chi-square test are said to be `overdispersed' with respect to
#'     the analytical uncertainties. The simplest way (`model-1') to
#'     deal with overdispersion is to inflate the standard error with
#'     a \eqn{\sqrt{MSWD}} premultiplier. In this case,
#'     \code{IsoplotR} returns two estimates of the standard error.
#'     To convert the second estimate to a confidence interval
#'     requires multiplication with the desired quantile of a
#'     t-distribution with the appropriate degrees of freedom.
#' @examples
#' attach(examples)
#' fit <- isochron(PbPb,plot=FALSE,exterr=FALSE)
#' err <- ci(x=fit$age[1],sx=fit$age[-1],oerr=5,df=fit$df)
#' message('age=',signif(fit$age[1],4),'Ma, ',
#'         '2se=',signif(err[1],2),'%, ',
#'         '2se(with dispersion)=',signif(err[2],2),'%')
#' @export
ci <- function(x,sx,oerr=3,df=NULL,absolute=FALSE){
    fact <- rep(ntfact(alpha()),length(sx))
    if (!is.null(df) && df>0){
        fact[-1] <- stats::qt(1-alpha()/2,df=df)
    }
    if (oerr>3 & absolute) oerr <- oerr-3
    if (oerr==1) out <- sx
    else if (oerr==2) out <- 2*sx
    else if (oerr==3) out <- fact*sx
    else if (oerr==4) out <- abs(100*sx/x)
    else if (oerr==5) out <- abs(200*sx/x)
    else if (oerr==6) out <- abs(100*fact*sx/x)
    else stop('Illegal oerr value')
    out
}

# formats table or vector of ages and errors
agerr <- function(x,...){ UseMethod("agerr",x) }
agerr.default <- function(x,oerr=1,sigdig=NA,...){
    out <- tst <- x
    tst[2] <- ci(x=x[1],sx=x[2],oerr=oerr)
    out <- roundit(age=tst[1],err=tst[2],sigdig=sigdig)
    out
}
agerr.matrix <- function(x,oerr=1,sigdig=NA,...){
    out <- tst <- x
    nc <- ncol(tst)
    i <- seq(from=2,to=nc,by=2)
    tst[,i] <- ci(x=x[,i-1],sx=x[,i],oerr=oerr)
    if (is.na(sigdig)){
        out <- tst
    } else {
        rounded <- roundit(age=tst[,i-1,drop=FALSE],
                           err=tst[,i,drop=FALSE],sigdig=sigdig)
        out[,i-1] <- rounded[,1:(nc/2)]
        out[,i] <- rounded[,(nc/2+1):nc]
        if (nc%%2==1) out[,nc] <- signif(tst[,nc],digits=sigdig)
    }
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
                    oerr=3,units=' Ma',prefix='age =',df=NULL){
    xerr <- ci(x,sx,oerr=oerr,df=df)
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
        if (identical(units,'%')) {
            werr <- ci(100*w,100*sw,oerr=oerr,absolute=TRUE)
            rounded <- roundit(100*w,werr,sigdig=sigdig,text=TRUE)
        } else {
            werr <- ci(w,sw,oerr=oerr)
            rounded <- roundit(w,werr,sigdig=sigdig,oerr=oerr,text=TRUE)
        }
    } else {
        w <- 0
        werr <- NA
        rounded <- c(w,werr)
    }
    lst <- list(p=prefix,a=rounded[1],b=rounded[2],u=units)
    if (oerr>3){
        out <- substitute(p~a*u%+-%b*'%',lst)
    } else {
        out <- substitute(p~a%+-%b*u,lst)
    }
    out
}
peaktit <- function(x,sx,p,sigdig=2,oerr=3,unit='Ma',prefix=NULL){
    xerr <- ci(x,sx,oerr=oerr)
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
get.ntit.default <- function(x,m=min(x,na.rm=TRUE),M=max(x,na.rm=TRUE),...){
    ns <- length(x)
    bad <- which(is.na(x) | x<m | x>M)
    nisnan <- length(bad)
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

# only used for measured disequilibrium
get.ci.ludwig <- function(par,x,type='joint',anchor=0,model=1,
                          exterr=FALSE,debug=FALSE){
    if (debug){ # test function
        if (TRUE){
            maxLL <- LL.ludwig(par,x=x,anchor=anchor,type=type,
                               model=model,exterr=exterr,LL48=FALSE,LL08=FALSE)
            pname <- 't'
            pval <- par['t']
        } else {
            maxLL <- LL.ludwig(par,x=x,anchor=anchor,type=type,
                               model=model,exterr=exterr,LL48=TRUE,LL08=FALSE)
            pname <- 'U48i'
            pval <- log(mclean(tt=exp(par['t']),d=x$d)$U48i)
        }
        pp <- seq(from=pval-.1,to=pval+.1,length.out=21)
        LL <- pp*0
        for (i in seq_along(pp)){
            LL[i] <- profile.LL.ludwig(pp[i],pname,x=x,maxLL=maxLL,
                                       type=type,anchor=anchor,
                                       model=model,exterr=exterr)
        }
        plot(pp,LL,type='b')
        lines(rep(pval,2),range(LL))
    }
    tt <- ifelse(anchor[1]==2,anchor[2],exp(par['t']))
    ci <- c(ll=NA,ul=NA)
    message('Computing profile likelihood of t')
    maxLL <- LL.ludwig(par,x=x,anchor=anchor,type=type,
                       model=model,exterr=exterr,LL48=FALSE,LL08=FALSE)
    ci['ll'] <- stats::uniroot(profile.LL.ludwig,interval=log(tt)-c(.1,0),
                               pname='t',x=x,maxLL=maxLL,type=type,
                               anchor=anchor,model=model,exterr=exterr,
                               extendInt="downX")$root
    ci['ul'] <- stats::uniroot(profile.LL.ludwig,interval=log(tt)+c(0,.1),
                               pname='t',x=x,maxLL=maxLL,type=type,anchor=anchor,
                               model=model,exterr=exterr,extendInt="upX")$root
    out <- data.frame('t'=ci)
    if (x$d$U48$option==2 && x$d$U48$sx>0){
        message('Computing profile likelihood of U48i')
        maxLL <- LL.ludwig(par,x=x,anchor=anchor,type=type,
                           model=model,exterr=exterr,LL48=TRUE,LL08=FALSE)
        U48i <- mclean(tt=tt,d=x$d)$U48i
        dLL <- profile.LL.ludwig(-Inf,pname='U48i',x=x,maxLL=maxLL,type=type,
                                 anchor=anchor,model=model,exterr=exterr)
        ci['ll'] <- ifelse(dLL<1.92,0,
                           stats::uniroot(profile.LL.ludwig,interval=log(U48i)-c(.1,0),
                                          pname='U48i',x=x,maxLL=maxLL,type=type,
                                          anchor=anchor,model=model,exterr=exterr,
                                          extendInt="downX")$root
                           )
        ci['ul'] <- stats::uniroot(profile.LL.ludwig,interval=log(U48i)+c(0,.1),
                                   pname='U48i',x=x,maxLL=maxLL,type=type,anchor=anchor,
                                   model=model,exterr=exterr,extendInt="upX")$root
        out[,'U48i'] <- ci
    }
    if (x$d$ThU$option==2 && x$d$ThU$sx>0){
        message('Computing profile likelihood of ThUi')
        maxLL <- LL.ludwig(par,x=x,anchor=anchor,type=type,
                           model=model,exterr=exterr,LL48=FALSE,LL08=TRUE)
        ThUi <- log(mclean(tt=tt,d=x$d)$ThUi)
        dLL <- profile.LL.ludwig(-Inf,pname='ThUi',x=x,maxLL=maxLL,type=type,
                                 anchor=anchor,model=model,exterr=exterr)
        ci['ll'] <- ifelse(dLL<1.92,0,
                           stats::uniroot(profile.LL.ludwig,interval=log(ThUi)-c(.1,0),
                                          pname='ThUi',x=x,maxLL=maxLL,type=type,
                                          anchor=anchor,model=model,exterr=exterr,
                                          extendInt="downX")$root
                           )
        ci['ul'] <- stats::uniroot(profile.LL.ludwig,interval=log(ThUi)+c(0,.1),
                                   pname='ThUi',x=x,maxLL=maxLL,type=type,anchor=anchor,
                                   model=model,exterr=exterr,extendInt="upX")$root
        out[,'U48i'] <- ci
    }
    exp(out)
}

profile.LL.ludwig <- function(thepar,pname,x,maxLL,type='joint',
                              anchor=0,model=1,exterr=FALSE,debug=FALSE){
    if (debug){
        browser()
    }
    X <- x
    if (pname=='t'){
        tt <- exp(thepar)
        anchor <- c(2,tt,0)
    } else if (pname=='U48i'){
        X$d$U48$x <- exp(thepar)
        X$d$U48$sx <- 0
        X$d$U48$option <- 1
    } else if (pname=='ThUi'){
        X$d$ThU$x <- exp(thepar)
        X$d$ThU$sx <- 0
        X$d$ThU$option <- 1
    } else {
        stop('No profile likelihood for parameter ',pname)
    }
    init <- init.ludwig(X,model=model,anchor=anchor,type=type)
    fit <- optim(init,fn=LL.ludwig,method='L-BFGS-B',
                 lower=init-2,upper=init+2,hessian=FALSE,x=X,
                 anchor=anchor,type=type,model=model,exterr=exterr)
    LL <- fit$value
    if (pname=='U48i'){
        pred <- mclean(tt=exp(fit$par['t']),d=X$d)
        LL <- LL - stats::dnorm(x=pred$U48,mean=x$d$U48$x,
                                sd=x$d$U48$sx,log=TRUE)
    } else if (pname=='ThUi'){
        pred <- mclean(tt=exp(fit$par['t']),d=X$d)
        LL <- LL - stats::dnorm(x=pred$ThU,mean=x$d$ThU$x,
                                sd=x$d$ThU$sx,log=TRUE)
    }
    abs(LL-maxLL)-1.92
}
