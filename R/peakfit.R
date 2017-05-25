#' Finite mixture modelling of geochronological datasets
#'
#' Implements the discrete mixture modelling algorithms of Galbraith
#' and Green (1993) and applies them to fission track and other
#' geochronological datasets.
#'
#' @param x either a \code{[2 x n]} matrix with measurements and their
#'     standard errors, or an object of class \code{fissiontracks},
#'     \code{UPb}, \code{ArAr}, \code{ReOs}, \code{SmNd}, \code{RbSr},
#'     \code{LuHf} or \code{UThHe}
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
#' @param ... optional arguments (not used)
#' @return a list with the following items: \describe{ \item{peaks}{a
#'     vector of peak locations} \item{props}{a vector of peak
#'     proportions} \item{peaks.err}{the standard errors of the peak
#'     locations} \item{props.err}{the standard errors of the peak
#'     proportions} \item{legend}{a vector of text expressions to be
#'     used in a figure legend} }
#' @references Galbraith, R.F. and Laslett, G.M., 1993. Statistical
#'     models for mixed fission track ages. Nuclear tracks and
#'     radiation measurements, 21(4), pp.459-470.
#' @examples
#' data(examples)
#' peakfit(examples$FT1,k=2)
#' @rdname peakfit
#' @export
peakfit <- function(x,...){ UseMethod("peakfit",x) }
#' @rdname peakfit
#' @export
peakfit.default <- function(x,k='auto',sigdig=2,log=TRUE,...){
    if (k<1) return(NULL)
    if (log) {
        x[,2] <- x[,2]/x[,1]
        x[,1] <- log(x[,1])
    }
    if (identical(k,'min')) {
        out <- min.age.model(x,sigdig=sigdig)
    } else if (identical(k,'auto')) {
        out <- normal.mixtures(x,k=BIC.fit(x,5),sigdig=sigdig,...)
    } else {
        out <- normal.mixtures(x,k,sigdig=sigdig,...)
    }
    if (log) {
        out$peaks <- exp(out$peaks)
        out$peaks.err <- out$peaks * out$peaks.err
    }
    out$legend <- peaks2legend(out,sigdig=sigdig,k=k)
    out
}
#' @rdname peakfit
#' @export
peakfit.fissiontracks <- function(x,k=1,exterr=TRUE,sigdig=2,log=TRUE,...){
    out <- NULL
    if (k == 0) return(out)
    if (identical(k,'auto')) k <- BIC.fit(x,5,log=log)
    if (x$format == 1 & !identical(k,'min')){
        out <- binomial.mixtures(x,k,exterr=exterr,...)
    }  else if (x$format == 3){
        out <- ages2peaks(x,k,log=log)
    } else {
        out <- peakfit.helper(x,k=k,sigdig=sigdig,log=log,
                              exterr=exterr,...)
    }
    out$legend <- peaks2legend(out,sigdig=sigdig,k=k)
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
                        sigdig=2,log=TRUE,...){
    peakfit.helper(x,k=k,type=type,cutoff.76=cutoff.76,
                   cutoff.disc=cutoff.disc,exterr=exterr,
                   sigdig=sigdig,log=log,...)
}
#' @param i2i
#'     `isochron to intercept': calculates the initial (aka `inherited',
#'     `excess', or `common') \eqn{^{40}}Ar/\eqn{^{36}}Ar,
#'     \eqn{^{87}}Sr/\eqn{^{86}}Sr, \eqn{^{143}}Nd/\eqn{^{144}}Nd,
#'     \eqn{^{187}}Os/\eqn{^{188}}Os or \eqn{^{176}}Hf/\eqn{^{177}}Hf
#'     ratio from an isochron fit. Setting \code{i2i} to \code{FALSE}
#'     uses the default values stored in \code{settings('iratio',...)}
#' @rdname peakfit
#' @export
peakfit.ArAr <- function(x,k=1,exterr=TRUE,sigdig=2,log=TRUE,i2i=FALSE,...){
    peakfit.helper(x,k=k,exterr=exterr,sigdig=sigdig,log=log,i2i=i2i,...)
}
#' @rdname peakfit
#' @export
peakfit.ReOs <- function(x,k=1,exterr=TRUE,sigdig=2,log=TRUE,i2i=TRUE,...){
    peakfit.helper(x,k=k,exterr=exterr,sigdig=sigdig,log=log,i2i=i2i,...)
}
#' @rdname peakfit
#' @export
peakfit.SmNd <- function(x,k=1,exterr=TRUE,sigdig=2,log=TRUE,i2i=TRUE,...){
    peakfit.helper(x,k=k,exterr=exterr,sigdig=sigdig,log=log,i2i=i2i,...)
}
#' @rdname peakfit
#' @export
peakfit.RbSr <- function(x,k=1,exterr=TRUE,sigdig=2,log=TRUE,i2i=TRUE,...){
    peakfit.helper(x,k=k,exterr=exterr,sigdig=sigdig,log=log,i2i=i2i,...)
}
#' @rdname peakfit
#' @export
peakfit.LuHf <- function(x,k=1,exterr=TRUE,sigdig=2,log=TRUE,i2i=TRUE,...){
    peakfit.helper(x,k=k,exterr=exterr,sigdig=sigdig,log=log,i2i=i2i,...)
}
#' @rdname peakfit
#' @export
peakfit.UThHe <- function(x,k=1,sigdig=2,log=TRUE,...){
    peakfit.helper(x,k=k,sigdig=sigdig,log=log,...)
}
peakfit.helper <- function(x,k=1,type=4,cutoff.76=1100,
                           cutoff.disc=c(-15,5),exterr=TRUE,sigdig=2,
                           log=TRUE,i2i=FALSE,...){
    if (k<1) return(NULL)
    if (identical(k,'auto')) k <- BIC.fit(x,5,log=log)
    fit <- ages2peaks(x,k=k,log=log,i2i=i2i)
    if (exterr){
        if (identical(k,'min')) numpeaks <- 1
        else numpeaks <- k
        for (i in 1:numpeaks){
            age.with.exterr <- add.exterr(x,fit$peaks[i],fit$peaks.err[i])
            fit$peaks.err[i] <- age.with.exterr[2]
        }
    }
    fit$legend <- peaks2legend(fit,sigdig=sigdig,k=k)
    fit    
}

ages2peaks <- function(x,k=1,type=4,cutoff.76=1100,
                       cutoff.disc=c(-15,5),log=TRUE,i2i=FALSE){
    if (hasClass(x,'UPb')){
        tt <- filter.UPb.ages(x,type,cutoff.76,
                              cutoff.disc,exterr=FALSE)
    } else if (hasClass(x,'ArAr')){
        tt <- ArAr.age(x,exterr=FALSE,i2i=i2i)
    } else if (hasClass(x,'UThHe')){
        tt <- UThHe.age(x)
    } else if (hasClass(x,'ReOs')){
        tt <- ReOs.age(x,exterr=FALSE,i2i=i2i)
    } else if (hasClass(x,'SmNd')){
        tt <- SmNd.age(x,exterr=FALSE,i2i=i2i)
    } else if (hasClass(x,'RbSr')){
        tt <- RbSr.age(x,exterr=FALSE,i2i=i2i)
    } else if (hasClass(x,'LuHf')){
        tt <- LuHf.age(x,exterr=FALSE,i2i=i2i)
    } else if (hasClass(x,'fissiontracks')){
        tt <- fissiontrack.age(x,exterr=FALSE)
    }
    peakfit.default(tt,k,log=log)
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

peaks2legend <- function(fit,sigdig=2,k=NULL){
    if (identical(k,'min')) return(minage2legend(fit,sigdig))
    out <- NULL
    for (i in 1:length(fit$peaks)){
        rounded.age <- roundit(fit$peaks[i],fit$peaks.err[i],sigdig=sigdig)
        rounded.prop <- roundit(fit$props[i],fit$props.err[i],sigdig=sigdig)
        line <- paste0('Peak ',i,': ',rounded.age[1],'\u00B1',
                       rounded.age[2],' (',100*rounded.prop[1],'\u00B1',
                       100*rounded.prop[2],'%)')
        out <- c(out,line)
    }
    out
}
minage2legend <- function(fit,sigdig=2){
    rounded.age <- roundit(fit$peaks,fit$peaks.err,sigdig=sigdig)
    paste0('Minimum: ',rounded.age[1],'+/-',rounded.age[2])
}

normal.mixtures <- function(x,k,sigdig=2,...){
    zu <- x[,1]
    su <- x[,2]
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
        fiu <- matrix(0,n,k)
        for (i in 1:k){
            fiu[,i] <- stats::dnorm(yu,betai[i]*xu,1)
        }
        piu <- matrix(0,n,k)
        for (u in 1:n){
            den <- sum(pii*fiu[u,])
            if (den>0) piu[u,] <- pii*fiu[u,]/den
        }
        fiu <- matrix(0,n,k)
        fu <- rep(0,n)
        for (i in 1:k){
            pii[i] <- mean(piu[,i])
            betai[i] <- sum(piu[,i]*xu*yu)/sum(piu[,i]*xu^2)
            fiu[,i] <- stats::dnorm(yu,betai[i]*xu,1)
            fu <- fu + pii[i] * fiu[,i]
        }
        newL <- sum(log(fu[fu>0]))
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
    peaks.err <- sqrt(diag(E)[k:(2*k-1)])
    props.err <- get.props.err(E)
    out <- list(L=L,peaks=betai,props=pii,
                peaks.err=peaks.err,props.err=props.err)
    out$legend <- peaks2legend(out,sigdig=sigdig,k=k)
    out
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
    props.err <- get.props.err(E)
    list(L=L,peaks=pe$peaks,props=pii,
         peaks.err=pe$peaks.err,props.err=props.err)
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

BIC.fit <- function(x,max.k,...){
    n <- numgrains(x)
    BIC <- Inf
    for (k in 1:max.k){
        fit <- peakfit(x,k,...)
        p <- 2*k-1
        newBIC <- -2*fit$L+p*log(n)
        if (newBIC<BIC) {
            BIC <- newBIC
        } else {
            k <- k-1
            break
        }
    }
    k
}

# Simple 3-parameter Normal model (Section 6.11 of Galbraith, 2005)
min.age.model <- function(zs,sigdig=2){
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
                newL <- get.min.age.L(pars,zs)
                if (newL < L) {
                    L <- newL
                    fit <- pars
                }
            }
        }
    }
    H <- stats::optimHess(fit,get.min.age.L,zs=zs)
    E <- solve(H)
    out <- list(L=L,peaks=fit[1],
                peaks.err=sqrt(E[1,1]),props.err=sqrt(E[3,3]))
    out
}

get.min.age.L <- function(pars,zs){
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
