#' @title
#' Finite mixture modelling of geochronological datasets
#' 
#' @description
#' Implements the discrete mixture modelling algorithms of Galbraith
#' and Laslett (1993) and applies them to fission track and other
#' geochronological datasets.
#' 
#' @details
#' Consider a dataset of \eqn{n} dates \eqn{\{t_1, t_2, ..., t_n\}}
#' with analytical uncertainties \eqn{\{s[t_1], s[t_2], ...,
#' s[t_n]\}}.  Define \eqn{z_i = \log(t_i)} and \eqn{s[z_i] =
#' s[t_i]/t_i}.  Suppose that these \eqn{n} values are derived from a
#' mixture of \eqn{k>2} populations with means
#' \eqn{\{\mu_1,...,\mu_k\}}. Such a \emph{discrete mixture} may be
#' mathematically described by \eqn{P(z_i|\mu,\omega) = \sum_{j=1}^k
#' \pi_j N(z_i | \mu_j, s[z_j]^2 )} where \eqn{\pi_j} is the
#' proportion of the population that belongs to the \eqn{j^{th}}
#' component, and \eqn{\pi_k=1-\sum_{j=1}^{k-1}\pi_j}. This equation
#' can be solved by the method of maximum likelihood (Galbraith and
#' Laslett, 1993).  \code{IsoplotR} implements the Bayes Information
#' Criterion (BIC) as a means of automatically choosing \eqn{k}. This
#' option should be used with caution, as the number of peaks steadily
#' rises with sample size (\eqn{n}).  If one is mainly interested in
#' the youngest age component, then it is more productive to use an
#' alternative parameterisation, in which all grains are assumed to
#' come from one of two components, whereby the first component is a
#' single discrete age peak (\eqn{\exp(m)}, say) and the second
#' component is a continuous distribution (as descibed by the
#' \code{\link{central}} age model), but truncated at this discrete
#' value. \code{IsoplotR} uses a normal likelihood function
#' (section 6.11 of Galbraith, 2005) for the minimum age model.
#' This may result in some inaccuracy for young and/or uranium-poor
#' fission track samples. 
#' 
#' @param x either an \code{[nx2]} matrix with measurements and their
#'     standard errors, or an object of class \code{fissiontracks},
#'     \code{UPb}, \code{PbPb}, \code{ThPb}, \code{ArAr}, \code{KCa},
#'     \code{ReOs}, \code{SmNd}, \code{RbSr}, \code{LuHf}, \code{ThU}
#'     or \code{UThHe}
#' @param k the number of discrete age components to be
#'     sought. Setting this parameter to \code{'auto'} automatically
#'     selects the optimal number of components (up to a maximum of 5)
#'     using the Bayes Information Criterion (BIC).
#' @param exterr propagate the external sources of uncertainty into
#'     the component age errors?
#' @param sigdig number of significant digits to be used for any
#'     legend in which the peak fitting results are to be displayed.
#' @param oerr indicates whether the analytical uncertainties of the
#'     output are reported in the plot legend as:
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
#' @param log take the logs of the data before applying the mixture
#'     model?
#' @param np number of parameters for the minimum age model. Must be
#'     either 3 or 4.
#' @param ... optional arguments (not used)
#' @seealso \code{\link{radialplot}}, \code{\link{central}}
#' @return Returns a list with the following items:
#'
#' \describe{
#'
#' \item{peaks}{a \code{2 x k} matrix with the following rows:
#'
#' \code{t}: the ages of the \code{k} peaks
#'
#' \code{s[t]}: the standard errors of \code{t}
#'
#' }
#'
#' \item{props}{a \code{2 x k} matrix with the following rows:
#'
#' \code{p}: the proportions of the \code{k} peaks
#'
#' \code{s[p]}: the standard errors of \code{p}
#'
#' }
#'
#' \item{L}{the log-likelihood of the fit}
#'
#' \item{legend}{a vector of text expressions to be used in a figure
#'     legend}
#'
#' }
#' @examples
#' attach(examples)
#' peakfit(FT1,k=2)
#'
#' peakfit(LudwigMixture,k='min')
#' @references
#' Galbraith, R.F. and Laslett, G.M., 1993. Statistical models for
#' mixed fission track ages. Nuclear Tracks and Radiation
#' Measurements, 21(4), pp.459-470.
#'
#' Galbraith, R.F. 2005, Statistics for fission track
#' analysis. Chapman and Hall/CRC, 229p.
#' 
#' @rdname peakfit
#' @export
peakfit <- function(x,...){ UseMethod("peakfit",x) }
#' @rdname peakfit
#' @export
peakfit.default <- function(x,k='auto',sigdig=2,oerr=3,log=TRUE,np=4,...){
    good <- !is.na(x[,1]+x[,2])
    X <- subset(x,subset=good)
    if (k<1) return(NULL)
    if (log) {
        X[,2] <- X[,2]/X[,1]
        X[,1] <- log(X[,1])
    }
    if (identical(k,'min')) {
        out <- min_age_model(X,np=np)
    } else if (identical(k,'auto')) {
        out <- normal.mixtures(X,k=BIC_fit(X,5),...)
    } else {
        out <- normal.mixtures(X,k,...)
    }
    if (log) {
        out$peaks['t',] <- exp(out$peaks['t',])
        out$peaks['s[t]',] <- out$peaks['t',] * out$peaks['s[t]',]
    }
    out$legend <- peaks2legend(out,k=k,sigdig=sigdig,oerr=oerr)
    out
}
#' @rdname peakfit
#' @export
peakfit.fissiontracks <- function(x,k=1,exterr=TRUE,sigdig=2,
                                  log=TRUE,oerr=3,np=4,...){
    out <- NULL
    if (k == 0) return(out)
    if (identical(k,'auto')) k <- BIC_fit(x,5,log=log)
    if (x$format == 1 & !identical(k,'min')){
        out <- binomial.mixtures(x,k,exterr=exterr,...)
    }  else if (x$format == 3){
        tt <- get.ages(x)
        out <- peakfit.default(tt,k=k,log=log,sigdig=sigdig,oerr=oerr)
    } else {
        out <- peakfit_helper(x,k=k,sigdig=sigdig,oerr=oerr,
                              log=log,exterr=exterr,...)
    }
    out
}
#' @param type scalar indicating whether to plot the
#'     \eqn{^{207}}Pb/\eqn{^{235}}U age (\code{type}=1), the
#'     \eqn{^{206}}Pb/\eqn{^{238}}U age (\code{type}=2), the
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb age (\code{type}=3), the
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb-\eqn{^{206}}Pb/\eqn{^{238}}U age
#'     (\code{type}=4), the concordia age (\code{type}=5), or the
#'     \eqn{^{208}}Pb/\eqn{^{232}}Th age (\code{type}=6).
#' @param cutoff.76 the age (in Ma) below which the
#'     \eqn{^{206}}Pb/\eqn{^{238}}U and above which the
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb age is used. This parameter is
#'     only used if \code{type=4}.
#' @param cutoff.disc discordance cutoff filter. This is an object of
#'     class \code{\link{discfilter}}.
#' 
#' @param common.Pb common lead correction:
#'
#' \code{0}:none
#'
#' \code{1}: use the Pb-composition stored in
#' 
#' \code{settings('iratio','Pb206Pb204')} (if \code{x} has class
#' \code{UPb} and \code{x$format<4});
#' 
#' \code{settings('iratio','Pb206Pb204')} and
#' \code{settings('iratio','Pb207Pb204')} (if \code{x} has class
#' \code{PbPb} or \code{x} has class \code{UPb} and
#' \code{3<x$format<7}); or
#'
#' \code{settings('iratio','Pb208Pb206')} and
#' \code{settings('iratio','Pb208Pb207')} (if \code{x} has class
#' \code{UPb} and \code{x$format=7} or \code{8}).
#'
#' \code{2}: use the isochron intercept as the initial Pb-composition
#'
#' \code{3}: use the Stacey-Kramers two-stage model to infer the
#' initial Pb-composition (only applicable if \code{x} has class
#' \code{UPb})
#' @rdname peakfit
#' @export
peakfit.UPb <- function(x,k=1,type=4,cutoff.76=1100,
                        cutoff.disc=discfilter(),common.Pb=0,
                        exterr=TRUE,sigdig=2,log=TRUE,oerr=3,np=4,...){
    peakfit_helper(x,k=k,type=type,cutoff.76=cutoff.76,
                   cutoff.disc=cutoff.disc,exterr=exterr,
                   sigdig=sigdig,log=log,oerr=oerr,
                   common.Pb=common.Pb,np=np,...)
}
#' @rdname peakfit
#' @export
peakfit.PbPb <- function(x,k=1,exterr=TRUE,sigdig=2,
                         log=TRUE,common.Pb=0,oerr=3,np=4,...){
    peakfit_helper(x,k=k,exterr=exterr,sigdig=sigdig,log=log,
                   common.Pb=common.Pb,oerr=oerr,np=np,...)
}
#' @param i2i `isochron to intercept': calculates the initial (aka
#'     `inherited', `excess', or `common')
#'     \eqn{^{40}}Ar/\eqn{^{36}}Ar, \eqn{^{40}}Ca/\eqn{^{44}}Ca,
#'     \eqn{^{207}}Pb/\eqn{^{204}}Pb, \eqn{^{87}}Sr/\eqn{^{86}}Sr,
#'     \eqn{^{143}}Nd/\eqn{^{144}}Nd, \eqn{^{187}}Os/\eqn{^{188}}Os,
#'     \eqn{^{230}}Th/\eqn{^{232}}Th, \eqn{^{176}}Hf/\eqn{^{177}}Hf or
#'     \eqn{^{204}}Pb/\eqn{^{208}}Pb ratio from an isochron
#'     fit. Setting \code{i2i} to \code{FALSE} uses the default values
#'     stored in \code{settings('iratio',...)}.
#' @rdname peakfit
#' @export
peakfit.ArAr <- function(x,k=1,exterr=TRUE,sigdig=2,
                         log=TRUE,i2i=FALSE,oerr=3,np=4,...){
    peakfit_helper(x,k=k,exterr=exterr,sigdig=sigdig,
                   log=log,i2i=i2i,oerr=oerr,np=np,...)
}
#' @rdname peakfit
#' @export
peakfit.ThPb <- function(x,k=1,exterr=TRUE,sigdig=2,
                         log=TRUE,i2i=FALSE,oerr=3,np=np,...){
    peakfit_helper(x,k=k,exterr=exterr,sigdig=sigdig,
                   log=log,i2i=i2i,oerr=oerr,np=4,...)
}
#' @rdname peakfit
#' @export
peakfit.KCa <- function(x,k=1,exterr=TRUE,sigdig=2,
                        log=TRUE,i2i=FALSE,oerr=3,np=4,...){
    peakfit_helper(x,k=k,exterr=exterr,sigdig=sigdig,
                   log=log,i2i=i2i,oerr=oerr,np=np,...)
}
#' @rdname peakfit
#' @export
peakfit.ReOs <- function(x,k=1,exterr=TRUE,sigdig=2,
                         log=TRUE,i2i=TRUE,oerr=3,np=4,...){
    peakfit_helper(x,k=k,exterr=exterr,sigdig=sigdig,
                   log=log,i2i=i2i,oerr=oerr,np=np,...)
}
#' @rdname peakfit
#' @export
peakfit.SmNd <- function(x,k=1,exterr=TRUE,sigdig=2,
                         log=TRUE,i2i=TRUE,oerr=3,np=4,...){
    peakfit_helper(x,k=k,exterr=exterr,sigdig=sigdig,
                   log=log,i2i=i2i,oerr=oerr,np=np,...)
}
#' @rdname peakfit
#' @export
peakfit.RbSr <- function(x,k=1,exterr=TRUE,sigdig=2,
                         log=TRUE,i2i=TRUE,oerr=3,np=4,...){
    peakfit_helper(x,k=k,exterr=exterr,sigdig=sigdig,
                   log=log,i2i=i2i,oerr=oerr,np=np,...)
}
#' @rdname peakfit
#' @export
peakfit.LuHf <- function(x,k=1,exterr=TRUE,sigdig=2,
                         log=TRUE,i2i=TRUE,oerr=3,np=4,...){
    peakfit_helper(x,k=k,exterr=exterr,sigdig=sigdig,
                   log=log,i2i=i2i,oerr=oerr,np=np,...)
}
#' @param Th0i initial \eqn{^{230}}Th correction.
#'
#' \code{0}: no correction
#'
#' \code{1}: project the data along an isochron fit
#'
#' \code{2}: if \code{x$format} is \code{1} or \code{2}, correct the
#' data using the measured present day \eqn{^{230}}Th/\eqn{^{238}}U,
#' \eqn{^{232}}Th/\eqn{^{238}}U and \eqn{^{234}}U/\eqn{^{238}}U
#' activity ratios in the detritus. If \code{x$format} is \code{3} or
#' \code{4}, correct the data using the measured
#' \eqn{^{238}}U/\eqn{^{232}}Th activity ratio of the whole rock, as
#' stored in \code{x} by the \code{read.data()} function.
#'
#' \code{3}: correct the data using an assumed initial
#' \eqn{^{230}}Th/\eqn{^{232}}Th-ratio for the detritus (only relevant
#' if \code{x$format} is \code{1} or \code{2}).
#' 
#' @rdname peakfit
#' @export
peakfit.ThU <- function(x,k=1,exterr=FALSE,sigdig=2,
                        log=TRUE,oerr=3,Th0i=0,np=4,...){
    peakfit_helper(x,k=k,exterr=exterr,sigdig=sigdig,
                   log=log,oerr=oerr,Th0i=Th0i,np=np,...)
}
#' @rdname peakfit
#' @export
peakfit.UThHe <- function(x,k=1,sigdig=2,log=TRUE,oerr=3,np=4,...){
    peakfit_helper(x,k=k,sigdig=sigdig,log=log,oerr=oerr,np=np,...)
}
peakfit_helper <- function(x,k=1,type=4,cutoff.76=1100,cutoff.disc=discfilter(),
                           exterr=TRUE,sigdig=2,log=TRUE,i2i=FALSE,
                           common.Pb=0,oerr=3,Th0i=0,np=4,...){
    if (k<1) return(NULL)
    if (identical(k,'auto')){
        k <- BIC_fit(x,5,log=log,type=type,cutoff.76=cutoff.76,i2i=i2i,
                     cutoff.disc=cutoff.disc,Th0i=Th0i,common.Pb=common.Pb)
    }
    tt <- get.ages(x,i2i=i2i,common.Pb=common.Pb,type=type,
                   cutoff.76=cutoff.76,cutoff.disc=cutoff.disc,Th0i=Th0i)
    fit <- peakfit.default(tt,k=k,np=np,log=log,oerr=oerr,sigdig=sigdig)
    if (exterr){
        if (identical(k,'min')) numpeaks <- 1
        else numpeaks <- k
        for (i in 1:numpeaks){
            age.with.exterr <- add.exterr(x,fit$peaks['t',i],
                                          fit$peaks['s[t]',i],type=type)
            fit$peaks['s[t]',i] <- age.with.exterr[2]
        }
    }
    fit$legend <- peaks2legend(fit,k=k,sigdig=sigdig,oerr=oerr)
    fit
}

get.peakfit.covmat <- function(k,pii,piu,aiu,biu){
    Au <- matrix(0,k-1,k-1)
    Bu <- matrix(0,k-1,k)
    Cu <- matrix(0,k,k)
    if (k>1){
        for (i in 1:(k-1)){
            for (j in 1:(k-1)){
                Au[i,j] <- sum((piu[,i]/pii[i]-piu[,k]/pii[k])*
                               (piu[,j]/pii[j]-piu[,k]/pii[k]))
            }
        }
        for (i in 1:(k-1)){
            for (j in 1:k){
                Bu[i,j] <- sum(piu[,j]*aiu[,j]*
                               (piu[,i]/pii[i]-piu[,k]/pii[k]-
                                (i==j)/pii[j]+(j==k)/pii[k])
                               )
            }
        }
    }
    for (i in 1:k){
        for (j in 1:k){
            Cu[i,j] <- sum(piu[,i]*piu[,j]*aiu[,i]*aiu[,j]-
                           (i==j)*biu[,i]*piu[,i])
        }
    }
    if (k>1) out <- blockinverse(Au,Bu,t(Bu),Cu,doall=TRUE)
    else out <- solve(Cu)
    out
}

get.props.err <- function(E){
    vars <- diag(E)
    k <- (nrow(E)+1)/2
    J <- -1 * matrix(1,1,k-1)
    if (k>1) {
        prop.k.err <- sqrt(J %*% E[1:(k-1),1:(k-1)] %*% t(J))
        out <- c(sqrt(vars[1:(k-1)]),prop.k.err)
    } else {
        out <- 0
    }
    out
}

peaks2legend <- function(fit,k=NULL,sigdig=2,oerr=3){
    if (identical(k,'min')){
        tit <- peaktit(x=fit$peaks[1],sx=fit$peaks[2],p=fit$props[1],
                       sigdig=sigdig,oerr=oerr,prefix=paste0('Minimum:'))
        out <- as.expression(tit)
    } else {
        out <- NULL
        for (i in 1:ncol(fit$peaks)){
            tit <- peaktit(x=fit$peaks[1,i],sx=fit$peaks[2,i],
                           p=fit$props[1,i],sigdig=sigdig,
                           oerr=oerr,prefix=paste0('Peak ',i,':'))
            out <- c(out,as.expression(tit))
        }
    }
    out
}

normal.mixtures <- function(x,k,...){
    good <- !is.na(x[,1]+x[,2])
    zu <- x[good,1]
    su <- x[good,2]
    xu <- 1/su
    yu <- zu/su
    n <- length(yu)
    if (k>1)
        betai <- seq(min(zu),max(zu),length.out=k)
    else
        betai <- stats::median(zu)
    pii <- rep(1,k)/k
    L <- -Inf
    for (j in 1:100){
        lpfiu <- matrix(0,n,k) # log(pii x fiu) taking logs enhances stability 
        for (i in 1:k){        # compared to Galbraith's original formulation
            lpfiu[,i] <- log(pii[i]) + stats::dnorm(yu,betai[i]*xu,1,log=TRUE)
        }
        piu <- matrix(0,n,k)
        for (i in 1:k){
            lden <- apply(lpfiu,2,'-',lpfiu[,i])
            piu[,i] <- 1/rowSums(exp(lden))
        }
        for (i in 1:k){
            pii[i] <- mean(piu[,i])
            betai[i] <- sum(piu[,i]*xu*yu)/sum(piu[,i]*xu^2)
            lpfiu[,i] <- log(pii[i]) +
                stats::dnorm(yu,betai[i]*xu,1,log=TRUE)
        }
        newL <- get.L.normal.mixture(lpfiu)
        if (((newL-L)/newL)^2 < 1e-20) break;
        L <- newL
    }
    aiu <- matrix(0,n,k)
    biu <- matrix(0,n,k)
    for (i in 1:k){
        aiu[,i] <- xu*(yu-betai[i]*xu)
        biu[,i] <- -(1-(yu-betai[i]*xu)^2)*xu^2
    }
    E <- get.peakfit.covmat(k,pii,piu,aiu,biu)
    out <- format.peaks(peaks=betai,peaks.err=sqrt(diag(E)[k:(2*k-1)]),
                        props=pii,props.err=get.props.err(E),df=n-2*k+1)
    out$L <- L
    out
}
# uses log-of-sums identity from Wikipedia
get.L.normal.mixture <- function(lpfiu){
    if (ncol(lpfiu)<2){
        fu <- lpfiu
    } else {
        sorted.lpfiu <- t(apply(lpfiu,1,sort,decreasing=TRUE))
        a0 <- subset(sorted.lpfiu,select=1)
        ai <- subset(sorted.lpfiu,select=-1)
        aia0 <- apply(ai,2,'-',a0)
        fu <- a0 + log(1 + rowSums(exp(aia0)))
    }
    sum(fu)
}

binomial.mixtures <- function(x,k,exterr=TRUE,...){
    yu <- x$x[,'Ns']
    mu <- x$x[,'Ns'] + x$x[,'Ni']
    NsNi <- (x$x[,'Ns']+0.5)/(x$x[,'Ni']+0.5)
    theta <- NsNi/(1+NsNi)
    thetai <- seq(min(theta),max(theta),length.out=k)
    pii <- rep(1,k)/k
    n <- length(yu)
    piu <- matrix(0,n,k)
    fiu <- matrix(0,n,k)
    aiu <- matrix(0,n,k)
    biu <- matrix(0,n,k)
    L <- -Inf
    # loop until convergence has been achieved
    for (j in 1:100){
        for (i in 1:k){
            fiu[,i] <- stats::dbinom(yu,mu,thetai[i])
        }
        for (u in 1:n){
            piu[u,] <- pii*fiu[u,]/sum(pii*fiu[u,])
        }
        fu <- rep(0,n)
        for (i in 1:k){
            pii[i] <- mean(piu[,i])
            thetai[i] <- sum(piu[,i]*yu)/sum(piu[,i]*mu)
            fu <- fu + pii[i] * fiu[,i]
        }
        newL <- sum(log(fu))
        if (((newL-L)/newL)^2 < 1e-20) break;
        L <- newL
    }
    for (i in 1:k){
        aiu[,i] <- yu-thetai[i]*mu
        biu[,i] <- (yu-thetai[i]*mu)^2 - thetai[i]*(1-thetai[i])*mu
    }
    E <- get.peakfit.covmat(k,pii,piu,aiu,biu)
    beta.var <- diag(E)[k:(2*k-1)]
    pe <- theta2age(x,thetai,beta.var,exterr)
    out <- format.peaks(peaks=pe$peaks,peaks.err=pe$peaks.err,
                        props=pii,props.err=get.props.err(E),df=n-2*k+1)
    out$L <- L
    out
}

format.peaks <- function(peaks,peaks.err,props,props.err,df){
    out <- list()
    k <- length(peaks)
    out$peaks <- matrix(0,2,k)
    colnames(out$peaks) <- 1:k
    rownames(out$peaks) <- c('t','s[t]')
    out$peaks['t',] <- peaks
    out$peaks['s[t]',] <- peaks.err
    out$props <- matrix(0,2,k)
    colnames(out$props) <- 1:k
    rownames(out$props) <- c('p','s[p]')
    out$props['p',] <- props
    out$props['s[p]',] <- props.err
    out
}

theta2age <- function(x,theta,beta.var,exterr=TRUE){
    rhoD <- x$rhoD
    zeta <- x$zeta
    if (!exterr) {
        rhoD[2] <- 0
        zeta[2] <- 0
    }
    k <- length(theta)
    peaks <- rep(0,k)
    peaks.err <- rep(0,k)
    for (i in 1:k){
        NsNi <- theta[i]/(1-theta[i])
        L8 <- lambda('U238')[1]
        peaks[i] <- log(1+0.5*L8*(zeta[1]/1e6)*rhoD[1]*(NsNi))/L8
        peaks.err[i] <- peaks[i]*sqrt(beta.var[i] +
                        (rhoD[2]/rhoD[1])^2 + (zeta[2]/zeta[1])^2)
    }
    list(peaks=peaks,peaks.err=peaks.err)
}

BIC_fit <- function(x,max.k,...){
    n <- length(x)
    BIC <- Inf
    tryCatch({
        for (k in 1:max.k){
            fit <- peakfit(x,k,...)
            p <- 2*k-1
            newBIC <- -2*fit$L+p*log(n)
            if (newBIC<BIC) {
                BIC <- newBIC
            } else {
                return(k-1)
            }
        }
        return(k)
    },error = function(e){
        return(k-1)
    }) 
}

min_age_model <- function(zs,np=4){
    # maps the parameters from -Inf/+Inf to model space
    mappar <- function(par,Mz){
        np <- length(par)
        gam <- par[1]
        prop <- exp(par[2])/(1+exp(par[2]))
        sig <- exp(par[3])
        if (np<4) mu <- gam
        else mu <- gam + (Mz-gam)*exp(par[4])/(1+exp(par[4]))
        c(gam,prop,sig,mu)
    }
    # computes the log-likelihood function using the transformed parameters
    LL <- function(par,zs,Mz){
        get.minage.L(pars=mappar(par,Mz),zs=zs)
    }
    mz <- min(zs[,1])
    Mz <- max(zs[,1])
    cfit <- continuous_mixture(zs[,1],zs[,2])
    init <- c(mz,0,log(cfit$sigma[1]),0)[1:np]
    ll <- c(mz,-20,init[3]-20,-20)[1:np]
    ul <- c(cfit$mu[1],20,init[3]+2,20)[1:np]
    fit <- stats::optim(init,LL,method='L-BFGS-B',zs=zs,Mz=Mz,lower=ll,upper=ul)
    if (fit$par[3]<(-10)) fit$par[2] <- 10 # sigma=0 => pi=1 (no sensitivity for pi)
    H <- tryCatch(stats::optimHess(fit$par,LL,zs=zs,Mz=Mz),
                  error=function(e){ return(NULL) })
    if (is.null(H) & np>3){
        return(min_age_model(zs=zs,np=3))
    }
    lE <- MASS::ginv(H)
    par <- mappar(fit$par,Mz)
    # propagate the uncertainties from -Inf/+Inf to model space
    J <- diag(np)
    J[2,2] <- exp(fit$par[2])/(1+exp(fit$par[2]))^2
    J[3,3] <- exp(fit$par[3])
    if (np==4){
        J[4,1] <- 1-exp(fit$par[4])/(1+exp(fit$par[4]))
        J[4,4] <- (Mz-par[1])*exp(fit$par[4])/(1+exp(fit$par[4]))^2
    }
    E <- J %*% lE %*% t(J)
    if (E[1,1]<0 & np==4){
        out <- min_age_model(zs=zs,np=3)
    } else {
        out <- list()
        out$L <- fit$value
        out$peaks <- matrix(0,2,1)
        rownames(out$peaks) <- c('t','s[t]')
        out$peaks['t',] <- par[1]
        out$peaks['s[t]',] <- if (E[1,1]<0) NA else sqrt(E[1,1])
        out$props <- matrix(0,2,1)
        rownames(out$props) <- c('p','s[p]')
        out$props['p',] <- par[2]
        out$props['s[p]',] <- if (E[2,2]<0) NA else sqrt(E[2,2])
        out$disp <- matrix(0,2,1)
        rownames(out$disp) <- c('d','s[d]')
        out$disp['d',] <- par[3]
        out$disp['s[d]',] <- if (E[3,3]<0) NA else sqrt(E[3,3])
        if (np==4){
            out$mu <- matrix(0,2,1)
            rownames(out$mu) <- c('t','s[t]')
            out$mu['t',] <- par[4]
            out$mu['s[t]',] <- if (E[4,4]<0) NA else sqrt(E[4,4])
        }
    }
    out
}

# Section 6.11 of Galbraith (2005)
get.minage.L <- function(pars,zs){
    z <- zs[,1]
    s <- zs[,2]
    gam <- pars[1]
    prop <- pars[2]
    sig <- pars[3]
    mu <- pars[4]
    AA  <- prop/sqrt(2*pi*s^2)
    BB <- -0.5*((z-gam)/s)^2
    CC <- (1-prop)/sqrt(2*pi*(sig^2+s^2))
    mu0 <- (mu/sig^2 + z/s^2)/(1/sig^2 + 1/s^2)
    s0 <- 1/sqrt(1/sig^2 + 1/s^2)
    DD <- 1-stats::pnorm((gam-mu0)/s0)
    EE <- 1-stats::pnorm((gam-mu)/sig)
    FF <- -0.5*((z-mu)^2)/(sig^2+s^2)
    fu <- AA*exp(BB) + CC*(DD/EE)*exp(FF)
    fu[fu<.Machine$double.xmin] <- .Machine$double.xmin
    fu[fu>.Machine$double.xmax] <- .Machine$double.xmax
    sum(-log(fu))
}
