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
#' mathematically described by:
#' \eqn{P(z_i|\mu,\omega) = \sum_{j=1}^k \pi_j N(z_i | \mu_j, s[z_j]^2 )}
#' where \eqn{\pi_j} is the proportion of the population that belongs
#' to the \eqn{j^{th}} component, and
#' \eqn{\pi_k=1-\sum_{j=1}^{k-1}\pi_j}. This equation can be solved by
#' the method of maximum likelihood (Galbraith and Laslett, 1993).
#' \code{IsoplotR} implements the Bayes Information Criterion (BIC) as
#' a means of automatically choosing \eqn{k}. This option should be
#' used with caution, as the number of peaks steadily rises with
#' sample size (\eqn{n}).  If one is mainly interested in the youngest
#' age component, then it is more productive to use an alternative
#' parameterisation, in which all grains are assumed to come from one
#' of two components, whereby the first component is a single discrete
#' age peak (\eqn{\exp(m)}, say) and the second component is a
#' continuous distribution (as descibed by the \code{\link{central}}
#' age model), but truncated at this discrete value (Van der Touw et
#' al., 1997).
#'
#' @param x either an \code{[n x 2]} matrix with measurements and
#'     their standard errors, or an object of class
#'     \code{fissiontracks}, \code{UPb}, \code{PbPb}, \code{ArAr},
#'     \code{KCa}, \code{ReOs}, \code{SmNd}, \code{RbSr}, \code{LuHf},
#'     \code{ThU} or \code{UThHe}
#' @param k the number of discrete age components to be
#'     sought. Setting this parameter to \code{'auto'} automatically
#'     selects the optimal number of components (up to a maximum of 5)
#'     using the Bayes Information Criterion (BIC).
#' @param exterr propagate the external sources of uncertainty into
#'     the component age errors?
#' @param sigdig number of significant digits to be used for any
#'     legend in which the peak fitting results are to be displayed.
#' @param log take the logs of the data before applying the mixture
#'     model?
#' @param alpha cutoff value for confidence intervals
#' @param ... optional arguments (not used)
#' @seealso \code{\link{radialplot}}, \code{\link{central}}
#' @return Returns a list with the following items:
#'
#' \describe{
#'
#' \item{peaks}{a \code{3 x k} matrix with the following rows:
#'
#' \code{t}: the ages of the \code{k} peaks
#'
#' \code{s[t]}: the estimated uncertainties of \code{t}
#'
#' \code{ci[t]}: the widths of approximate \eqn{100(1-\alpha)\%}
#' confidence intervals for \code{t}}
#'
#' \item{props}{a \code{2 x k} matrix with the following rows:
#'
#' \code{p}: the proportions of the \code{k} peaks
#'
#' \code{s[p]}: the estimated uncertainties (standard errors) of
#' \code{p}}
#'
#' \item{L}{the log-likelihood of the fit}
#'
#' \item{legend}{a vector of text expressions to be used in a figure
#'     legend}
#'
#' }
#' @references
#' Galbraith, R.F. and Laslett, G.M., 1993. Statistical models for
#' mixed fission track ages. Nuclear Tracks and Radiation
#' Measurements, 21(4), pp.459-470.
#'
#' van der Touw, J., Galbraith, R., and Laslett, G. A logistic
#' truncated normal mixture model for overdispersed binomial
#' data. Journal of Statistical Computation and Simulation,
#' 59(4):349-373, 1997.
#' @examples
#' data(examples)
#' peakfit(examples$FT1,k=2)
#'
#' peakfit(examples$LudwigMixture,k='min')
#' @rdname peakfit
#' @export
peakfit <- function(x,...){ UseMethod("peakfit",x) }
#' @rdname peakfit
#' @export
peakfit.default <- function(x,k='auto',sigdig=2,log=TRUE,alpha=0.05,...){
    good <- !is.na(x[,1]+x[,2])
    x <- subset(x,subset=good)
    if (k<1) return(NULL)
    if (log) {
        x[,2] <- x[,2]/x[,1]
        x[,1] <- log(x[,1])
    }
    if (identical(k,'min')) {
        out <- min_age_model(x,sigdig=sigdig,alpha=alpha)
    } else if (identical(k,'auto')) {
        out <- normal.mixtures(x,k=BIC_fit(x,5),sigdig=sigdig,alpha=alpha,...)
    } else {
        out <- normal.mixtures(x,k,sigdig=sigdig,alpha=alpha,...)
    }
    if (log) {
        out$peaks['t',] <- exp(out$peaks['t',])
        out$peaks['s[t]',] <- out$peaks['t',] * out$peaks['s[t]',]
        out$peaks['ci[t]',] <- out$peaks['t',] * out$peaks['ci[t]',]
    }
    out$legend <- peaks2legend(out,sigdig=sigdig,k=k)
    out
}
#' @rdname peakfit
#' @export
peakfit.fissiontracks <- function(x,k=1,exterr=TRUE,sigdig=2,
                                  log=TRUE,alpha=0.05,...){
    out <- NULL
    if (k == 0) return(out)
    if (identical(k,'auto')) k <- BIC_fit(x,5,log=log)
    if (x$format == 1 & !identical(k,'min')){
        out <- binomial.mixtures(x,k,exterr=exterr,alpha=alpha,...)
    }  else if (x$format == 3){
        tt <- get.ages(x)
        out <- peakfit.default(tt,k=k,log=log,alpha=alpha)
    } else {
        out <- peakfit_helper(x,k=k,sigdig=sigdig,log=log,
                              exterr=exterr,alpha=alpha,...)
    }
    out$legend <- peaks2legend(out,sigdig=sigdig,k=k)
    out
}
#'
#' @param type scalar valueindicating whether to plot the
#'     \eqn{^{207}}Pb/\eqn{^{235}}U age (\code{type}=1), the
#'     \eqn{^{206}}Pb/\eqn{^{238}}U age (\code{type}=2), the
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb age (\code{type}=3), the
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb-\eqn{^{206}}Pb/\eqn{^{238}}U age
#'     (\code{type}=4), or the (Wetherill) concordia age
#'     (\code{type}=5)
#' 
#' @param cutoff.76 the age (in Ma) below which the
#'     \eqn{^{206}}Pb/\eqn{^{238}}U and above which the
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb age is used. This parameter is
#'     only used if \code{type=4}.
#' 
#' @param cutoff.disc two element vector with the maximum and minimum
#'     percentage discordance allowed between the
#'     \eqn{^{207}}Pb/\eqn{^{235}}U and \eqn{^{206}}Pb/\eqn{^{238}}U
#'     age (if \eqn{^{206}}Pb/\eqn{^{238}}U < \code{cutoff.76}) or
#'     between the \eqn{^{206}}Pb/\eqn{^{238}}U and
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb age (if
#'     \eqn{^{206}}Pb/\eqn{^{238}}U > \code{cutoff.76}).  Set
#'     \code{cutoff.disc=NA} if you do not want to use this filter.
#' 
#' @param common.Pb apply a common lead correction using one of three
#'     methods:
#'
#' \code{1}: use the Stacey-Kramers two-stage model to infer the initial
#' Pb-composition
#'
#' \code{2}: use the isochron intercept as the initial Pb-composition
#'
#' \code{3}: use the Pb-composition stored in
#' \code{settings('iratio','Pb207Pb206')} (if \code{x$format}<4) or
#' \code{settings('iratio','Pb206Pb204')} and
#' \code{settings('iratio','Pb207Pb204')} (if \code{x$format}>3)
#' 
#' @rdname peakfit
#' @export
peakfit.UPb <- function(x,k=1,type=4,cutoff.76=1100,cutoff.disc=c(-15,5),
                        common.Pb=0,exterr=TRUE,sigdig=2,log=TRUE,alpha=0.05,...){
    peakfit_helper(x,k=k,type=type,cutoff.76=cutoff.76,cutoff.disc=cutoff.disc,
                   exterr=exterr,sigdig=sigdig,log=log,alpha=alpha,common.Pb=common.Pb,...)
}
#' @rdname peakfit
#' @export
peakfit.PbPb <- function(x,k=1,exterr=TRUE,sigdig=2,
                         log=TRUE,common.Pb=0,alpha=0.05,...){
    peakfit_helper(x,k=k,exterr=exterr,sigdig=sigdig,
                   log=log,common.Pb=common.Pb,alpha=alpha,...)
}
#'
#' @param i2i `isochron to intercept': calculates the initial (aka
#'     `inherited', `excess', or `common')
#'     \eqn{^{40}}Ar/\eqn{^{36}}Ar, \eqn{^{40}}Ca/\eqn{^{44}}Ca,
#'     \eqn{^{87}}Sr/\eqn{^{86}}Sr, \eqn{^{143}}Nd/\eqn{^{144}}Nd,
#'     \eqn{^{187}}Os/\eqn{^{188}}Os or \eqn{^{176}}Hf/\eqn{^{177}}Hf
#'     ratio from an isochron fit. Setting \code{i2i} to \code{FALSE}
#'     uses the default values stored in
#'     \code{settings('iratio',...)}. When applied to data of class
#'     \code{ThU}, setting \code{i2i} to \code{TRUE} applies a
#'     detrital Th-correction.
#' 
#' @rdname peakfit
#' @export
peakfit.ArAr <- function(x,k=1,exterr=TRUE,sigdig=2,
                         log=TRUE,i2i=FALSE,alpha=0.05,...){
    peakfit_helper(x,k=k,exterr=exterr,sigdig=sigdig,
                   log=log,i2i=i2i,alpha=alpha,...)
}
#' @rdname peakfit
#' @export
peakfit.KCa <- function(x,k=1,exterr=TRUE,sigdig=2,
                         log=TRUE,i2i=FALSE,alpha=0.05,...){
    peakfit_helper(x,k=k,exterr=exterr,sigdig=sigdig,
                   log=log,i2i=i2i,alpha=alpha,...)
}
#' @rdname peakfit
#' @export
peakfit.ReOs <- function(x,k=1,exterr=TRUE,sigdig=2,
                         log=TRUE,i2i=TRUE,alpha=0.05,...){
    peakfit_helper(x,k=k,exterr=exterr,sigdig=sigdig,
                   log=log,i2i=i2i,alpha=alpha,...)
}
#' @rdname peakfit
#' @export
peakfit.SmNd <- function(x,k=1,exterr=TRUE,sigdig=2,
                         log=TRUE,i2i=TRUE,alpha=0.05,...){
    peakfit_helper(x,k=k,exterr=exterr,sigdig=sigdig,
                   log=log,i2i=i2i,alpha=alpha,...)
}
#' @rdname peakfit
#' @export
peakfit.RbSr <- function(x,k=1,exterr=TRUE,sigdig=2,
                         log=TRUE,i2i=TRUE,alpha=0.05,...){
    peakfit_helper(x,k=k,exterr=exterr,sigdig=sigdig,
                   log=log,i2i=i2i,alpha=alpha,...)
}
#' @rdname peakfit
#' @export
peakfit.LuHf <- function(x,k=1,exterr=TRUE,sigdig=2,
                         log=TRUE,i2i=TRUE,alpha=0.05,...){
    peakfit_helper(x,k=k,exterr=exterr,sigdig=sigdig,
                   log=log,i2i=i2i,alpha=alpha,...)
}
#' @param detritus detrital \eqn{^{230}}Th correction (only applicable
#'     when \code{x$format == 1} or \code{2}.
#'
#' \code{0}: no correction
#'
#' \code{1}: project the data along an isochron fit
#'
#' \code{2}: correct the data using an assumed initial
#' \eqn{^{230}}Th/\eqn{^{232}}Th-ratio for the detritus.
#'
#' \code{3}: correct the data using the measured present day
#' \eqn{^{230}}Th/\eqn{^{238}}U, \eqn{^{232}}Th/\eqn{^{238}}U and
#' \eqn{^{234}}U/\eqn{^{238}}U-ratios in the detritus.
#'
#' @rdname peakfit
#' @export
peakfit.ThU <- function(x,k=1,exterr=FALSE,sigdig=2, log=TRUE,
                        i2i=TRUE,alpha=0.05,detritus=0,...){
    peakfit_helper(x,k=k,exterr=exterr,sigdig=sigdig,log=log,i2i=i2i,
                   alpha=alpha,detritus=detritus,...)
}
#' @rdname peakfit
#' @export
peakfit.UThHe <- function(x,k=1,sigdig=2,log=TRUE,alpha=0.05,...){
    peakfit_helper(x,k=k,sigdig=sigdig,log=log,alpha=alpha,...)
}
peakfit_helper <- function(x,k=1,type=4,cutoff.76=1100,cutoff.disc=c(-15,5),
                           exterr=TRUE,sigdig=2,log=TRUE,i2i=FALSE,
                           common.Pb=0,alpha=0.05,detritus=0,...){
    if (k<1) return(NULL)
    if (identical(k,'auto')){
        k <- BIC_fit(x,5,log=log,type=type,cutoff.76=cutoff.76,
                     cutoff.disc=cutoff.disc,detritus=detritus,commonPb=common.Pb)
    }
    tt <- get.ages(x,i2i=i2i,common.Pb=common.Pb,type=type,cutoff.76=cutoff.76,
                   cutoff.disc=cutoff.disc,detritus=detritus)
    fit <- peakfit.default(tt,k=k,sigdig=sigdig,log=log,alpha=alpha)
    if (exterr){
        if (identical(k,'min')) numpeaks <- 1
        else numpeaks <- k
        for (i in 1:numpeaks){
            age.with.exterr <- add.exterr(x,fit$peaks['t',i],fit$peaks['s[t]',i],type=type)
            fit$peaks['s[t]',i] <- age.with.exterr[2]
            fit$peaks['ci[t]',i] <- nfact(alpha)*fit$peaks['s[t]',i]
        }
    }
    fit$legend <- peaks2legend(fit,sigdig=sigdig,k=k)
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

peaks2legend <- function(fit,sigdig=2,k=NULL){
    if (identical(k,'min')) return(min_age_to_legend(fit,sigdig))
    out <- NULL
    for (i in 1:ncol(fit$peaks)){
        rounded.age <- roundit(fit$peaks[1,i],fit$peaks[2:3,i],sigdig=sigdig)
        line <- paste0('Peak ',i,': ',rounded.age[1],' \u00B1 ',
                       rounded.age[2],' | ',rounded.age[3])
        if (k>1){
            rounded.prop <- roundit(100*fit$props[1,i],100*fit$props[2,i],sigdig=sigdig)
            line <- paste0(line,' (',rounded.prop[1],'%)')
        }
        out <- c(out,line)
    }
    out
}
min_age_to_legend <- function(fit,sigdig=2){
    rounded.age <- roundit(fit$peaks[1],fit$peaks[2:3],sigdig=sigdig)
    paste0('Minimum: ',rounded.age[1],'\u00B1',rounded.age[2],' | ',rounded.age[3])
}

normal.mixtures <- function(x,k,sigdig=2,alpha=0.05,...){
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
        for (i in 1:k){        # compared to Galbraiths orignal formulation
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
    out <- format.peaks(peaks=betai,
                        peaks.err=sqrt(diag(E)[k:(2*k-1)]),
                        props=pii,
                        props.err=get.props.err(E),
                        df=n-2*k+1,alpha=alpha)
    out$L <- L
    out$legend <- peaks2legend(out,sigdig=sigdig,k=k)
    out
}
# uses sum-of-logs identity from Wikipedia
get.L.normal.mixture <- function(lpfiu){
    if (ncol(lpfiu)<2){
        fu <- lpfiu
    } else {
        sorted.lpfiu <- t(apply(lpfiu,1,sort,decreasing=TRUE))
        a0 <- subset(sorted.lpfiu,select=1)
        ai <- subset(sorted.lpfiu,select=-1)
        aia0 <- apply(ai,2,'-',a0)
        fu <- a0 + log(1 + sum(exp(aia0)))
    }
    sum(fu)
}

binomial.mixtures <- function(x,k,exterr=TRUE,alpha=0.05,...){
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
    out <- format.peaks(peaks=pe$peaks,
                        peaks.err=pe$peaks.err,
                        props=pii,
                        props.err=get.props.err(E),
                        df=n-2*k+1,alpha=alpha)
    out$L <- L
    out
}

format.peaks <- function(peaks,peaks.err,props,props.err,df,alpha=0.05){
    out <- list()
    k <- length(peaks)
    out$peaks <- matrix(0,3,k)
    colnames(out$peaks) <- 1:k
    rownames(out$peaks) <- c('t','s[t]','ci[t]')
    out$peaks['t',] <- peaks
    out$peaks['s[t]',] <- peaks.err
    out$peaks['ci[t]',] <- nfact(alpha)*out$peaks['s[t]',]
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

BIC_fit <- function(x,max.k,type=4,cutoff.76=1100,cutoff.disc=c(-15,5),
                    exterr=TRUE,detritus=0,common.Pb=0,...){
    n <- length(x)
    BIC <- Inf
    tryCatch({
        for (k in 1:max.k){
            fit <- peakfit(x,k,type=type,cutoff.76=cutoff.76,
                           cutoff.disc=cutoff.disc,exterr=exterr,
                           detritus=detritus,common.Pb=common.Pb,...)
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

# Simple 3-parameter Normal model (Section 6.11 of Galbraith, 2005)
min_age_model <- function(zs,sigdig=2,alpha=0.05){
    z <- zs[,1]
    mu <- seq(min(z),max(z),length.out=100)
    sigma <- seq(stats::sd(z)/10,2*stats::sd(z),length.out=10)
    prop <- seq(0,1,length.out=20)
    L <- Inf
    # grid search!
    for (mui in mu){ # mu
        for (sigmai in sigma){ # sigma
            for (propi in prop){ # pi
                pars <- c(mui,sigmai,propi)
                newL <- get.minage.L(pars,zs)
                if (newL < L) {
                    L <- newL
                    fit <- pars
                }
            }
        }
    }
    H <- stats::optimHess(fit,get.minage.L,zs=zs)
    E <- solve(H)
    out <- list()
    out$L <- L
    out$peaks <- matrix(0,3,1)
    rownames(out$peaks) <- c('t','s[t]','ci[t]')
    df <- length(z)-3
    out$peaks['t',] <- fit[1]
    if (E[1,1]<0) out$peaks['s[t]',] <- NA
    else out$peaks['s[t]',] <- sqrt(E[1,1])
    out$peaks['ci[t]',] <- nfact(alpha)*out$peaks['s[t]',]
    out
}

get.minage.L <- function(pars,zs){
    z <- zs[,1]
    s <- zs[,2]
    mu <- pars[1]
    sigma <- pars[2]
    prop <- pars[3]
    AA  <- prop/sqrt(2*pi*s^2)
    BB <- -0.5*((z-mu)/s)^2
    CC <- (1-prop)/sqrt(s*pi*(sigma^2+s^2))
    mu0 <- (mu/sigma^2 + z/s^2)/(1/sigma^2 + 1/s^2)
    s0 <- 1/sqrt(1/sigma^2 + 1/s^2)
    DD <- 2*(1-stats::pnorm((mu-mu0)/s0))
    EE <- -0.5*((z-mu)^2)/(sigma^2+s^2)
    fu <- AA*exp(BB) + CC*DD*exp(EE)
    -sum(log(fu))
}
