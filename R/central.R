#' Calculate U-Th-He and fission track central ages and compositions
#'
#' Computes the geometric mean composition of a continuous mixture of
#' fission track or U-Th-He data and returns the corresponding age and
#' fitting parameters.
#'
#' @details
#' The central age assumes that the observed age distribution is the
#' combination of two sources of scatter: analytical uncertainty and
#' true geological dispersion.
#' \enumerate{
#' \item For fission track data, the analytical uncertainty is assumed
#' to obey Poisson counting statistics and the geological dispersion
#' is assumed to follow a lognormal distribution.
#' \item For U-Th-He data, the U-Th-(Sm)-He compositions and
#' uncertainties are assumed to follow a logistic normal distribution.
#' \item For all other data types, both the analytical uncertainties
#' and the true ages are assumed to follow lognormal distributions.
#' }
#' The difference between the central age and the weighted mean age is
#' usually small unless the data are imprecise and/or strongly
#' overdispersed.
#'
#' @param x an object of class \code{UThHe} or \code{fissiontracks},
#'     OR a 2-column matrix with (strictly positive) values and
#'     uncertainties
#' @param alpha cutoff value for confidence intervals
#' @param ... optional arguments
#' @return
#' If \code{x} has class \code{UThHe}, returns a list containing the
#'     following items:
#'
#' \describe{
#'
#' \item{uvw}{(if the input data table contains Sm) or \strong{uv} (if
#' it does not): the mean log[U/He], log[Th/He] (, and log[Sm/He])
#' composition.}
#'
#' \item{covmat}{the covariance matrix of \code{uvw} or \code{uv}.}
#'
#' \item{mswd}{the reduced Chi-square statistic of data concordance,
#' i.e. \eqn{mswd=SS/df}, where \eqn{SS} is the sum of squares of the
#' log[U/He]-log[Th/He] compositions.}
#'
#' \item{model}{the fitting model.}
#'
#' \item{df}{the degrees of freedom (\eqn{2n-2}) of the fit (only
#' reported if \code{model=1}).}
#'
#' \item{p.value}{the p-value of a Chi-square test with \code{df}
#' degrees of freedom (only reported if \code{model=1}.)}
#'
#' }
#'
#' \item{age}{a three- or four-element vector with:
#'
#' \code{t}: the central age.
#'
#' \code{s[t]}: the standard error of \code{t}.
#'
#' \code{ci[t]}: the width of a \eqn{100(1-\alpha)\%} confidence
#' interval for \code{t}.
#'
#' \code{disp[t]}: the studentised \eqn{100(1-\alpha)\%} confidence
#' interval enhanced by a factor of \eqn{\sqrt{mswd}} (only reported
#' if \code{model=1}).
#'
#' }
#'
#' \item{w}{the geological overdispersion term. If \code{model=3},
#' this is a three-element vector with the standard deviation of the
#' (assumedly) Normal dispersion and the lower and upper half-widths
#' of its \eqn{100(1-\alpha)\%} confidence interval. \code{w=0} if
#' code{model<3}.}
#'
#' OR, otherwise:
#'
#' \describe{
#'
#' \item{age}{a three-element vector with:
#'
#' \code{t}: the central age.
#'
#' \code{s[t]}: the standard error of \code{t}.
#'
#' \code{ci[t]}: the width of a \eqn{100(1-\alpha)\%} confidence
#' interval for \code{t}.}
#'
#' \item{disp}{a three-element vector with the overdispersion
#' (standard deviation) of the excess scatter, and the upper and lower
#' half-widths of its \eqn{100(1-\alpha)\%} confidence interval.}
#'
#' \item{mswd}{the reduced Chi-square statistic of data concordance,
#' i.e. \eqn{mswd=X^2/df}, where \eqn{X^2} is a Chi-square statistic
#' of the EDM data or ages}
#'
#' \item{df}{the degrees of freedom (\eqn{n-2})}
#'
#' \item{p.value}{the p-value of a Chi-square test with \code{df}
#' degrees of freedom}
#' }
#'
#' @seealso \code{\link{weightedmean}}, \code{\link{radialplot}},
#'     \code{\link{helioplot}}
#' @examples
#' data(examples)
#' print(central(examples$UThHe)$age)
#'
#' @references Galbraith, R.F. and Laslett, G.M., 1993. Statistical
#'     models for mixed fission track ages. Nuclear Tracks and
#'     Radiation Measurements, 21(4), pp.459-470.
#'
#' Vermeesch, P., 2008. Three new ways to calculate average (U-Th)/He
#'     ages. Chemical Geology, 249(3), pp.339-347.
#'
#' @rdname central
#' @export
central <- function(x,...){ UseMethod("central",x) }
#' @rdname central
#' @export
central.default <- function(x,alpha=0.05,...){
    good <- !is.na(rowSums(x))
    zu <- log(x[good,1])
    su <- x[good,2]/x[good,1]
    fit <- continuous_mixture(zu,su)
    tt <- exp(fit$mu[1])
    st <- tt*fit$mu[2]
    Chi2 <- sum((zu/su)^2,na.rm=TRUE)-(sum(zu/su^2,na.rm=TRUE)^2)/
        sum(1/su^2,na.rm=TRUE)
    out <- list()
    # remove two d.o.f. for mu and sigma
    out$df <- length(zu)-2
    # add back one d.o.f. for the homogeneity test
    out$mswd <- Chi2/(out$df+1)
    out$p.value <- 1-stats::pchisq(Chi2,out$df+1)
    out$age <- c(tt,st,nfact(alpha)*st)
    out$disp <- c(fit$sigma,
                  profile_LL_weightedmean_disp(fit,zu,su,alpha))
    names(out$age) <- c('t','s[t]','ci[t]')
    names(out$disp) <- c('s','ll','ul')
    out
}
#' @param model choose one of the following statistical models:
#'
#' \code{1}: weighted mean. This model assumes that the scatter
#' between the data points is solely caused by the analytical
#' uncertainty. If the assumption is correct, then the MSWD value
#' should be approximately equal to one. There are three strategies to
#' deal with the case where MSWD>1. The first of these is to assume
#' that the analytical uncertainties have been underestimated by a
#' factor \eqn{\sqrt{MSWD}}.
#'
#' \code{2}: unweighted mean. A second way to deal with over- or
#' underdispersed datasets is to simply ignore the analytical
#' uncertainties.
#'
#' \code{3}: weighted mean with overdispersion: instead of attributing
#' any overdispersion (MSWD > 1) to underestimated analytical
#' uncertainties (model 1), one could also attribute it to the
#' presence of geological uncertainty, which manifests itself as an
#' added (co)variance term.
#'
#' @rdname central
#' @export
central.UThHe <- function(x,alpha=0.05,model=1,...){
    out <- list()
    ns <- nrow(x)
    doSm <- doSm(x)
    fit <- UThHe_logratio_mean(x,model=model,w=0)
    mswd <- mswd_UThHe(x,fit,doSm=doSm)
    f <- tfact(alpha,mswd$df)
    if (model==1){
        out <- c(fit,mswd)
        out$age['disp[t]'] <-
            uvw2age(out,doSm=doSm(x),fact=f*sqrt(mswd$mswd))[2]
    } else if (model==2){
        out <- fit
    } else {
        f <- nfact(alpha)
        w <- get.UThHe.w(x,fit)
        out <- UThHe_logratio_mean(x,model=model,w=w)
        out$w <- c(w,profile_LL_central_disp_UThHe(fit=out,x=x,alpha=alpha))
        names(out$w) <- c('s','ll','ul')
    }
    out$age['ci[t]'] <- f*out$age['s[t]']
    out
}
#' @param mineral setting this parameter to either \code{apatite} or
#'     \code{zircon} changes the default efficiency factor, initial
#'     fission track length and density to preset values (only affects
#'     results if \code{x$format=2})
#' @rdname central
#' @export
central.fissiontracks <- function(x,mineral=NA,alpha=0.05,...){
    out <- list()
    if (x$format<2){
        L8 <- lambda('U238')[1]
        sigma <- 0.15 # convenient starting value
        Nsj <- x$x[,'Ns']
        Ns <- sum(Nsj)
        Nij <- x$x[,'Ni']
        Ni <- sum(Nij)
        num <- (Nsj*Ni-Nij*Ns)^2
        den <- Nsj+Nij
        Chi2 <- sum(num/den)/(Ns*Ni)
        mj <- Nsj+Nij
        pj <- Nsj/mj
        theta <- Ns/sum(mj)
        for (i in 1:30){ # page 49 of Galbraith (2005)
            wj <- mj/(theta*(1-theta)+(mj-1)*(theta*(1-theta)*sigma)^2)
            sigma <- sigma * sqrt(sum((wj*(pj-theta))^2)/sum(wj))
            theta <- sum(wj*pj)/sum(wj)
        }
        tt <- log(1+0.5*L8*(x$zeta[1]/1e6)*x$rhoD[1]*theta/(1-theta))/L8
        st <- tt * sqrt( 1/(sum(wj)*(theta*(1-theta))^2) +
                         (x$rhoD[2]/x$rhoD[1])^2 +
                         (x$zeta[2]/x$zeta[1])^2 )
        mu <- log(theta/(1-theta))
        # remove two d.o.f. for mu and sigma
        out$df <- length(Nsj)-2
        # add back one d.o.f. for homogeneity test
        out$mswd <- Chi2/(out$df+1)
        out$p.value <- 1-stats::pchisq(Chi2,out$df+1)
        out$age <- c(tt,st,stats::qt(1-alpha/2,out$df)*st)
        out$disp <- c(sigma,profile_LL_central_disp_FT(
                                mu=mu,sigma=sigma,y=Nsj,m=mj,alpha=alpha))
        names(out$age) <- c('t','s[t]','ci[t]')
        names(out$disp) <- c('s','ll','ul')
    } else if (x$format>1){
        tst <- age(x,exterr=FALSE,mineral=mineral)
        out <- central.default(tst,alpha=alpha)
    }
    out
}

UThHe_logratio_mean <- function(x,model=1,w=0,fact=1){
    out <- average_uvw(x,model=model,w=w)
    out$model <- model
    out$w <- w
    out$age <- rep(NA,3)
    names(out$age) <- c('t','s[t]','ci[t]')
    out$age[c('t','s[t]')] <- uvw2age(out,doSm(x))
    out
}

uvw2age <- function(fit,doSm,fact=1){
    if (doSm){
        cc <- uvw2UThHe(fit$uvw,fact*fit$covmat)
        out <- get.UThHe.age(cc['U'],cc['sU'],
                             cc['Th'],cc['sTh'],
                             cc['He'],cc['sHe'],
                             cc['Sm'],cc['sSm'])
    } else {
        cc <- uv2UThHe(fit$uvw,fact*fit$covmat)
        out <- get.UThHe.age(cc['U'],cc['sU'],
                             cc['Th'],cc['sTh'],
                             cc['He'],cc['sHe'])
    }
    out
}

average_uvw <- function(x,model=1,w=0){
    out <- list()
    doSm <- doSm(x)
    if (doSm(x)){
        nms <- c('u','v','w')
        logratios <- flat.uvw.table(x,w=w)
        fit <- wtdmean3D(logratios)
        uvw <- logratios[,c(1,3,5)]
    } else {
        nms <- c('u','v')
        logratios <- flat.uv.table(x,w=w)
        fit <- wtdmean2D(logratios)
        uvw <- logratios[,c(1,3)]
    }
    if (model==2){
        out$uvw <- apply(uvw,2,mean)
        out$covmat <- stats::cov(uvw)/(nrow(uvw)-1)
    } else {
        out$uvw <- fit$x
        out$covmat <- fit$cov
    }
    names(out$uvw) <- nms
    colnames(out$covmat) <- nms
    rownames(out$covmat) <- nms
    out
}

get.UThHe.w <- function(x,fit){
    stats::optimize(LL.uvw,interval=c(0,100),
                    UVW=fit$uvw,x=x,doSm=doSm(x),
                    maximum=TRUE)$maximum
}

LL.uvw <- function(w,UVW,x,doSm=TRUE,LL=TRUE){
    out <- 0
    for (i in 1:length(x)){
        if (doSm){
            uvwc <- UThHe2uvw.covmat(x,i,w=w)
            X <- UVW-uvwc$uvw
        } else {
            uvwc <- UThHe2uv.covmat(x,i,w=w)
            X <- UVW-uvwc$uv
        }
        E <- uvwc$covmat
        SS <- X %*% solve(E) %*% t(X)
        if (LL) {
            out <- out - ( determinant(E,logarithm=TRUE)$modulus + SS )/2
        } else {
            out <- out + SS
        }
    }
    out
}
mswd_UThHe <- function(x,fit,doSm=FALSE){
    out <- list()
    if (doSm) nd <- 3
    else nd <- 2
    SS <- LL.uvw(w=fit$w,fit$uvw[1:nd],x,doSm=doSm,LL=FALSE)
    out$df <- nd*(length(x)-1)
    out$mswd <- as.numeric(SS/out$df)
    out$p.value <- 1-stats::pchisq(SS,out$df)
    out
}

continuous_mixture <- function(zu,su){
    sigma <- stats::sd(zu)
    for (i in 1:30){ # page 100 of Galbraith (2005)
        wu <- 1/(sigma^2+su^2)
        mu <- sum(wu*zu,na.rm=TRUE)/sum(wu,na.rm=TRUE)
        fit <- stats::optimize(eq.6.9,c(0,10*stats::sd(zu)),mu=mu,zu=zu,su=su)
        sigma <- fit$minimum
    }
    smu <- 1/sqrt(sum(wu,na.rm=TRUE))
    out <- list()
    out$mu <- c(mu,smu)
    out$sigma <- sigma
    out
}

eq.6.9 <- function(sigma,mu,zu,su){
    wu <- 1/(sigma^2+su^2)
    (1-sum((wu*(zu-mu))^2,na.rm=TRUE)/sum(wu,na.rm=TRUE))^2
}
