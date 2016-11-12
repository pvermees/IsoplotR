#' Finite mixture modelling of geochronological datasets
#'
#' Implements the discrete mixture modelling algorithms of Galbraith
#' and Green (1993) and applies them to fission track and other
#' geochronological datasets.
#'
#' @param x either a \code{[2 x n]} matrix with measurements and their
#'     standard errors, or an object of class \code{fissiontracks}
#' @param k the number of discrete age components to be
#'     sought. Setting this parameter to \code{'auto'} automatically
#'     selects the optimal number of components (up to a maximum of 5)
#'     using the Bayes Information Criterion (BIC).
#' @param exterr propagate the external sources of uncertainty into
#'     the component age errors?
#' @param sigdig number of significant digits to be used for any legend
#' in which the peak fitting results are to be displayed.
#' @param ... optional arguments (not used)
#' @return a list with the following items:
#' \describe{
#' \item{peaks}{a vector of peak locations}
#' \item{props}{a vector of peak proportions}
#' \item{peaks.err}{the standard errors of the peak locations}
#' \item{props.err}{the standard errors of the peak proportions}
#' \item{legend}{a vector of text expressions to be used in a figure legend}
#' }
#' @references
#' Galbraith, R.F. and Laslett, G.M., 1993. Statistical models for
#' mixed fission track ages. Nuclear tracks and radiation
#' measurements, 21(4), pp.459-470.
#' @examples
#' data(examples)
#' peakfit(examples$FT1,k=2)
#' @rdname peakfit
#' @export
peakfit <- function(x,...){ UseMethod("peakfit",x) }
#' @rdname peakfit
#' @export
peakfit.default <- function(x,k=1,sigdig=2,...){
    if (k<1) return(NULL)
    zu <- x[,1]
    su <- x[,2]
    xu <- 1/su
    yu <- zu/su
    n <- length(yu)
    betai <- seq(min(zu),max(zu),length.out=k)
    pii <- rep(1,k)/k
    piu <- matrix(0,n,k)
    fiu <- matrix(0,n,k)
    aiu <- matrix(0,n,k)
    biu <- matrix(0,n,k)
    L <- -Inf
    for (j in 1:100){
        for (i in 1:k){
            fiu[,i] <- stats::dnorm(yu,betai[i]*xu,1)
        }
        for (u in 1:n){
            piu[u,] <- pii*fiu[u,]/sum(pii*fiu[u,])
        }
        fu <- rep(0,n)
        for (i in 1:k){
            pii[i] <- mean(piu[,i])
            betai[i] <- sum(piu[,i]*xu*yu)/sum(piu[,i]*xu^2)
            fu <- fu + pii[i] * fiu[,i]
        }
        newL <- sum(log(fu))
        if (((newL-L)/newL)^2 < 1e-20) break;
        L <- newL
    }
    for (i in 1:k){
        aiu[,i] <- xu*(yu-betai[i]*xu)
        biu[,i] <- -(1-(yu-betai[i]*xu)^2)*xu^2
    }
    E <- get.peakfit.covmat(k,pii,piu,aiu,biu)
    peaks.err <- sqrt(diag(E)[k:(2*k-1)])
    props.err <- get.props.err(E)
    out <- list(peaks=betai,props=pii,
                peaks.err=peaks.err,props.err=props.err)
    out$legend <- peaks2legend(out,sigdig=sigdig)
    out
}
#' @rdname peakfit
#' @export
peakfit.fissiontracks <- function(x,k=1,exterr=TRUE,sigdig=2,...){
    if (k == 0) return(NULL)
    if (x$format == 1){
        yu <- x$x[,'Ns']
        mu <- x$x[,'Ns'] + x$x[,'Ni']
        theta <- (x$x[,'Ns']/x$x[,'Ni'])/(1+x$x[,'Ns']/x$x[,'Ni'])
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
        theta.var <- diag(E)[k:(2*k-1)]
        peaks <- rep(0,k)
        peaks.err <- rep(0,k)
        rhoD <- x$rhoD
        zeta <- x$zeta
        if (!exterr) {
            rhoD[2] <- 0
            zeta[2] <- 0
        }
        for (i in 1:k){
            NsNi <- thetai[i]/(1-thetai[i])
            relErrNsNi <- theta.var[i]/thetai[i]^2
            L8 <- lambda('U238')[1]
            peaks[i] <- log(1+0.5*L8*(zeta[1]/1e6)*rhoD[1]*(NsNi))/L8
            peaks.err[i] <- peaks[i]*sqrt(relErrNsNi +
                            (rhoD[2]/rhoD[1])^2 + (zeta[2]/zeta[1])^2)
        }
        props.err <- get.props.err(E)
        out <- list(peaks=peaks,props=pii,
                    peaks.err=peaks.err,props.err=props.err)
    } else {
        out <- ages2peaks(x,k)
    }
    out$legend <- peaks2legend(out,sigdig=sigdig)
    out
}
#' @param type scalar indicating whether to plot the
#'     \eqn{^{207}}Pb/\eqn{^{235}}U age (\code{type}=1), the
#'     \eqn{^{206}}Pb/\eqn{^{238}}U age (\code{type}=2), the
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb age (\code{type}=3), the
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb-\eqn{^{206}}Pb/\eqn{^{238}}U age
#'     (\code{type}=4), or the (Wetherill) concordia age (\code{type}=5)
#' @param cutoff.76 the age (in Ma) below which the
#'     \eqn{^{206}}Pb/\eqn{^{238}}U and above which the
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb age is used. This parameter is
#'     only used if \code{type=4}.
#' @param cutoff.disc two element vector with the maximum and minimum
#'     percentage discordance allowed between the
#'     \eqn{^{207}}Pb/\eqn{^{235}}U and \eqn{^{206}}Pb/\eqn{^{238}}U
#'     age (if \eqn{^{206}}Pb/\eqn{^{238}}U < \code{cutoff.76}) or
#'     between the \eqn{^{206}}Pb/\eqn{^{238}}U and
#'     \eqn{^{207}}Pb/\eqn{^{206}}Pb age (if
#'     \eqn{^{206}}Pb/\eqn{^{238}}U > \code{cutoff.76}).  Set
#'     \code{cutoff.disc=NA} if you do not want to use this filter.
#' @rdname peakfit
#' @export
peakfit.UPb <- function(x,k=1,type=4,cutoff.76=1100,
                        cutoff.disc=c(-15,5),exterr=TRUE,
                        sigdig=2,...){
    if (k<1) return(NULL)
    fit <- ages2peaks(x,k=k,type=type,cutoff.76=cutoff.76,
                      cutoff.disc=cutoff.disc)
    if (exterr){
        for (i in 1:k){
            R <- get.ratios.UPb(tt=fit$peaks[i],st=fit$peaks.err[i],
                                exterr=TRUE,as.UPb=TRUE)
            age.with.exterr <- filter.UPb.ages(R,type=type,cutoff.76=cutoff.76,
                                               cutoff.disc=cutoff.disc,exterr=TRUE)
            fit$peaks.err[i] <- age.with.exterr[2]
        }
    }
    fit$legend <- peaks2legend(fit,sigdig=sigdig)
    fit
}
#' @rdname peakfit
#' @export
peakfit.ArAr <- function(x,k=1,exterr=TRUE,sigdig=2,...){
    if (k<1) return(NULL)
    fit <- ages2peaks(x,k=k)
    if (exterr){
        for (i in 1:k){
            R <- get.ArAr.ratio(fit$peaks[i],fit$peaks.err[i],
                                x$J[1],0,exterr=FALSE)
            age.with.exterr <- get.ArAr.age(R[1],R[2],x$J[1],x$J[2],exterr=TRUE)
            fit$peaks.err[i] <- age.with.exterr[2]
        }
    }
    fit$legend <- peaks2legend(fit,sigdig=sigdig)
    fit
}
#' @rdname peakfit
#' @export
peakfit.UThHe <- function(x,k=1,sigdig=2,...){
    if (k<1) return(NULL)
    fit <- ages2peaks(x,k=k)
    fit$legend <- peaks2legend(fit,sigdig=sigdig)
    fit
}

ages2peaks <- function(x,k=1,type=4,cutoff.76=1100,
                       cutoff.disc=c(-15,5),log=TRUE){
    if (hasClass(x,'UPb')){
        tt <- filter.UPb.ages(x,type,cutoff.76,
                              cutoff.disc,exterr=FALSE)
    } else if (hasClass(x,'ArAr')){
        tt <- ArAr.age(x,exterr=FALSE)
    } else if (hasClass(x,'UThHe')){
        tt <- UThHe.age(x)
    } else if (hasClass(x,'fissiontracks')){
        tt <- fissiontrack.age(x,exterr=FALSE)
    }
    if (log) zs <- cbind(log(tt[,1]),tt[,2]/tt[,1])
    else zs <- tt
    out <- peakfit.default(zs,k)
    if (log) {
        out$peaks <- exp(out$peaks)
        out$peaks.err <- out$peaks * out$peaks.err
    }
    out
}

get.peakfit.covmat <- function(k,pii,piu,aiu,biu){
    Au <- matrix(0,k-1,k-1)
    Bu <- matrix(0,k-1,k)
    Cu <- matrix(0,k,k)
    hess <- matrix(0,2*k-1,2*k-1)
    if (k>1){
        for (i in 1:(k-1)){
            for (j in 1:(k-1)){
                Au[i,j] <- sum((piu[,i]/pii[i]-piu[,k]/pii[k])*
                               (piu[,j]/pii[j]-piu[,k]/pii[k]))
            }
        }
        hess[1:(k-1),1:(k-1)] <- Au
        for (i in 1:(k-1)){
            for (j in 1:k){
                Bu[i,j] <- sum(piu[,j]*aiu[,j]*
                               (piu[,i]/pii[i]-piu[,k]/pii[k]-
                                (i==j)/pii[j]+(j==k)/pii[k])
                               )
            }
        }
        hess[1:(k-1),k:(2*k-1)] <- Bu
        hess[k:(2*k-1),1:(k-1)] <- t(Bu)
    }
    for (i in 1:k){
        for (j in 1:k){
            Cu[i,j] <- sum(piu[,i]*piu[,j]*aiu[,i]*aiu[,j]-
                           (i==j)*biu[,i]*piu[,i])
        }
    }
    hess[k:(2*k-1),k:(2*k-1)] <- Cu
    solve(hess)
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

peaks2legend <- function(fit,sigdig=2){
    out <- NULL
    for (i in 1:length(fit$peaks)){
        rounded.age <- roundit(fit$peaks[i],fit$peaks.err[i],sigdig=sigdig)
        rounded.prop <- roundit(fit$props[i],fit$props.err[i],sigdig=sigdig)
        line <- paste0('Peak ',i,': ',rounded.age$x,'+/-',
                       rounded.age$err,' (',rounded.prop$x,'+/-',
                       rounded.prop$err,'%)')
        out <- c(out,line)
    }
    out
}
