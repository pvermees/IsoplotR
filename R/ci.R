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
ci <- function(x=0,sx,oerr=3,df=NULL,absolute=FALSE){
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

errlabelswapper <- function(errnames){
    if (any(grepl('s\\[',errnames))){
        out <- gsub('s\\[','err\\[',errnames)
    } else {
        out <- errnames
    }
    out
}

#' Formats table or vector of ages and errors
#' @param x an IsoplotR data object
#' @noRd
agerr <- function(x,...){ UseMethod("agerr",x) }
#' @param oerr output error (1:4 for 1 or 2 sigma, absolute or relative)
#' @param sigdig number of significant digits
#' @noRd
agerr.default <- function(x,oerr=1,sigdig=NA,...){
    out <- tst <- x
    tst[2] <- ci(x=x[1],sx=x[2],oerr=oerr)
    out <- roundit(age=tst[1],err=tst[2],sigdig=sigdig)
    names(out) <- names(tst) <- errlabelswapper(names(x))
    out
}
#' @noRd
agerr.matrix <- function(x,oerr=1,sigdig=NA,...){
    out <- tst <- x
    nc <- ncol(x)
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
    colnames(out) <- errlabelswapper(colnames(x))
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
peaktit <- function(x,sx,p=NULL,sigdig=2,oerr=3,unit='Ma',prefix=NULL){
    xerr <- ci(x,sx,oerr=oerr)
    rounded.x <- roundit(x,xerr,sigdig=sigdig,oerr=oerr,text=TRUE)
    if (is.null(p)){
        lst <- list(p=prefix,a=rounded.x[1],b=rounded.x[2],u=unit)
        if (oerr>3){
            out <- substitute(p~a~u%+-%b,lst)
        } else {
            out <- substitute(p~a%+-%b~u,lst)
        }
    } else {
        rounded.p <- signif(100*p,sigdig)
        lst <- list(p=prefix,a=rounded.x[1],b=rounded.x[2],c=rounded.p[1],u=unit)
        if (oerr>3){
            out <- substitute(p~a~u%+-%b*'% (prop='*c*'%)',lst)
        } else {
            out <- substitute(p~a%+-%b~u~'(prop='*c*'%)',lst)
        }
    }
    out
}
bayestit <- function(x,XL,n=NULL,ntit=paste0('(n=',n,')'),
                     sigdig=2,oerr=3,units=' Ma',prefix='age ='){
    if (oerr%in%c(1,4)){
        lq <- -1
        uq <- 1
    } else if (oerr%in%c(2,5)){
        lq <- -2
        uq <- 2
    } else if (oerr%in%c(3,6)){
        lq <- stats::qnorm(alpha())
        uq <- stats::qnorm(1-alpha())
    }
    cdf <- cumsum(XL[,'L'])/sum(XL[,'L'])
    increasing <- which(diff(cdf)>1e-20)
    ll <- stats::spline(x=cdf[increasing],y=XL[increasing,'x'],
                        xout=stats::pnorm(lq),method='hyman')
    ul <- stats::spline(x=cdf[increasing],y=XL[increasing,'x'],
                        xout=stats::pnorm(uq),method='hyman')
    le <- (x-ll$y)
    ue <- (ul$y-x)
    if (oerr>3) {
        le <- 100*le/x
        ue <- 100*ue/x
    }
    rounded <- roundit(x,c(ue,le),sigdig=sigdig,oerr=oerr,text=TRUE)
    lst <- list(p=prefix,a=rounded[1],b=rounded[2],c=rounded[3],u=units,n=ntit)
    if (oerr>3){
        out <- substitute(p~a*u+b/-c*'%'~n,lst)
    } else {
        out <- substitute(p~a+b/-c*u~n,lst)
    }
    out
}

#' Get a sample size label for the plot title
#' @param x an IsoplotR data object or a data vector
#' @noRd
get_ntit <- function(x,...){
    UseMethod("get_ntit",x)
}
#' @param m minimum cutoff. Values below this are removed from the
#'     total sample size
#' @param M maximum cutoff. Valus above this are removed from the
#'     total sample size
#' @noRd
get_ntit.default <- function(x,m=min(x,na.rm=TRUE),M=max(x,na.rm=TRUE),...){
    ns <- length(x)
    bad <- which(is.na(x) | x<m | x>M)
    nisnan <- length(bad)
    out <- '(n='
    if (nisnan>0) out <- paste0(out,ns-nisnan,'/')
    paste0(out,ns,')')
}
#' @noRd
get_ntit.fissiontracks <- function(x,...){
    if (x$format<2){
        out <- get_ntit.default(x$x[,'Ns'])
    } else {
        out <- get_ntit.default(x$Ns)
    }
    out    
}
