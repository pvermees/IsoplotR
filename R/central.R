#' @title
#' Fits random effects models to overdispersed datasets
#'
#' @description
#' Computes the logratio mean composition of a continuous mixture of
#' fission track or U-Th-He data and returns the corresponding age and
#' fitting parameters. Only propagates the systematic uncertainty
#' associated with decay constants and calibration factors after
#' computing the weighted mean isotopic composition. Does not propagate
#' the uncertainty of any initial daughter correction, because this is
#' neither a purely random or purely systematic uncertainty.
#'
#' @details
#' The central age assumes that the observed age distribution is the
#' combination of two sources of scatter: analytical uncertainty and
#' true geological dispersion.
#' \enumerate{
#' \item For fission track data, the analytical uncertainty is assumed
#' to obey Binomial counting statistics and the geological dispersion
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
#' The uncertainty budget of the central age does not include the
#' uncertainty of the initial daughter correction (if any), for the
#' same reasons as discussed under the \code{\link{weightedmean}}
#' function.
#' 
#' @param x an object of class \code{UThHe} or \code{fissiontracks},
#'     OR a 2-column matrix with (strictly positive) values and
#'     uncertainties
#' @param ... optional arguments
#' @return If \code{x} has class \code{UThHe} and \code{compositional}
#'     is \code{TRUE}, returns a list containing the following items:
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
#' \item{age}{a two- or three-element vector with:\cr
#' \code{t}: the 'barycentric' age, i.e. the age corresponding to \code{uvw}.\cr
#' \code{s[t]}: the standard error of \code{t}.\cr
#' \code{disp[t]}: the standard error of \code{t} enhanced by a
#' factor of \eqn{\sqrt{mswd}} (only reported if \code{model=1}). }
#'
#' \item{w}{the geological overdispersion term. If \code{model=3},
#' this is a two-element vector with the standard deviation of the
#' (assumedly) Normal dispersion and its standard error. \code{w=0} if
#' \code{model<3}.}
#'
#' }
#'
#' OR, otherwise:
#'
#' \describe{
#'
#' \item{age}{a two-element vector with the central age and its
#' standard error.}
#'
#' \item{disp}{a two-element vector with the overdispersion (standard
#' deviation) of the excess scatter, and its standard error.}
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
#' 
#' @examples
#' attach(examples)
#' print(central(UThHe)$age)
#'
#' @references
#'
#' Galbraith, R.F. and Laslett, G.M., 1993. Statistical models for
#' mixed fission track ages. Nuclear Tracks and Radiation
#' Measurements, 21(4), pp.459-470.
#'
#' Vermeesch, P., 2008. Three new ways to calculate average (U-Th)/He
#'     ages. Chemical Geology, 249(3), pp.339-347.
#'
#' @rdname central
#' @export
central <- function(x,...){ UseMethod("central",x) }
#' @rdname central
#' @export
central.default <- function(x,...){
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
    out <- append(out,getMSWD(Chi2,out$df+1))
    out$age <- c(tt,st)
    out$disp <- fit$sigma
    names(out$age) <- c('t','s[t]')
    names(out$disp) <- c('w','s[w]')
    out
}
#' @param compositional logical. If \code{TRUE}, calculates the
#'     'barycentric' U-Th-He, age, i.e. the age corresponding to the
#'     weighted mean logratio composition.
#' @param model only relevant if \code{compositional} is
#'     \code{TRUE}. If the scatter between the data points is solely
#'     caused by the analytical uncertainty, then the MSWD value
#'     should be approximately equal to one. There are three
#'     strategies to deal with the case where MSWD>1.choose one of the
#'     following statistical models:
#'
#' \code{1}: assume that the analytical uncertainties have been
#' underestimated by a factor \eqn{\sqrt{MSWD}}.
#'
#' \code{2}: ignore the analytical
#' uncertainties.
#'
#' \code{3}: attribute any excess dispersion to the presence of
#' geological uncertainty, which manifests itself as an added
#' (co)variance term.
#'
#' @rdname central
#' @export
central.UThHe <- function(x,compositional=FALSE,model=1,...){
    if (compositional){
        ns <- nrow(x)
        doSm <- doSm(x)
        fit <- UThHe_logratio_mean(x,model=model,w=0)
        mswd <- mswd_UThHe(x,fit,doSm=doSm)
        mswd$model <- model # to fulfil requirements of inflate function
        if (inflate(mswd)){
            fit$age['disp[t]'] <- uvw2age(fit,doSm=doSm,fact=mswd$mswd)[2]
        }
        if (model==1){
            out <- c(fit,mswd)
        } else if (model==2){
            out <- fit
        } else {
            init <- log(sqrt(mswd$mswd)*fit$age['s[t]']/fit$age['t'])
            lw <- stats::optimise(LL.uvw,interval=init+c(-5,5),
                                  UVW=fit$uvw,x=x,doSm=doSm(x),
                                  maximum=TRUE)$maximum
            w <- exp(lw)
            H <- stats::optimHess(lw,LL.uvw,UVW=fit$uvw,x=x,doSm=doSm(x))
            sw <- w*solve(-H)
            out <- UThHe_logratio_mean(x,model=model,w=w)
            out$disp <- c(w,sw)
            names(out$disp) <- c('w','s[w]')
        }
    } else {
        out <- central.default(age(x))
    }
    out
}

#' @param exterr include the zeta or decay constant uncertainty into
#'     the error propagation for the central age?
#' @rdname central
#' @export
central.fissiontracks <- function(x,exterr=FALSE,...){
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
        st <- tt * sqrt(1/(sum(wj)*(theta*(1-theta))^2) +
                        (x$rhoD[2]/x$rhoD[1])^2 + (x$zeta[2]/x$zeta[1])^2)
        mu <- log(theta/(1-theta))
        # remove two d.o.f. for mu and sigma
        out$df <- length(Nsj)-2
        # add back one d.o.f. for homogeneity test
        out$mswd <- Chi2/(out$df+1)
        out$p.value <- 1-stats::pchisq(Chi2,out$df+1)
        out$age <- c(tt,st)
        kappa <- log(sigma)
        H <- stats::optimHess(par=c(mu,kappa),fn=LL.FT,y=Nsj,m=mj)
        out$disp <- c(sigma,sqrt(solve(-H)[2,2])*sigma)
        names(out$age) <- c('t','s[t]')
        names(out$disp) <- c('w','s[w]')
    } else if (x$format>1){
        tst <- age(x,exterr=FALSE)
        out <- central.default(tst)
    }
    if (exterr){
        out$age[1:2] <- add.exterr(x,tt=out$age[1],st=out$age[2])
    }
    out
}

UThHe_logratio_mean <- function(x,model=1,w=0,fact=1){
    out <- average_uvw(x,model=model,w=w)
    out$model <- model
    out$disp <- w
    out$age <- uvw2age(out,doSm(x),fact=fact)
    names(out$age) <- c('t','s[t]')
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
        uvw <- logratios[,c(1,3,5),drop=FALSE]
    } else {
        nms <- c('u','v')
        logratios <- flat.uv.table(x,w=w)
        fit <- wtdmean2D(logratios)
        uvw <- logratios[,c(1,3),drop=FALSE]
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

LL.uvw <- function(lw,UVW,x,doSm=TRUE,LL=TRUE){
    w <- exp(lw)
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
    as.numeric(out)
}
mswd_UThHe <- function(x,fit,doSm=FALSE){
    out <- list()
    if (doSm) nd <- 3
    else nd <- 2
    SS <- LL.uvw(lw=log(fit$disp),fit$uvw[1:nd],x,doSm=doSm,LL=FALSE)
    out$df <- nd*(length(x)-1)
    out$mswd <- as.numeric(SS/out$df)
    out$p.value <- 1-stats::pchisq(SS,out$df)
    out
}

continuous_mixture <- function(zu,su){
    if (length(unique(zu))>1){
        sigma <- stats::sd(zu)
        for (i in 1:30){ # page 100 of Galbraith (2005)
            wu <- 1/(sigma^2+su^2)
            mu <- sum(wu*zu,na.rm=TRUE)/sum(wu,na.rm=TRUE)
            kappa <- log(sigma)
            fit <- stats::optimise(f=eq.6.9,interval=kappa+c(-3,3),
                                   mu=mu,zu=zu,su=su)
            sigma <- exp(fit$minimum)
        }
    } else{
        sigma <- 0
        wu <- 1/su^2
        mu <- zu[1]
    }
    out <- list()
    smu <- 1/sqrt(sum(wu,na.rm=TRUE))
    ssigma <- 1/sqrt(2*(sigma^2)*sum(wu^2,na.rm=TRUE))
    out$mu <- c(mu,smu)
    out$sigma <- c(sigma,ssigma)
    out
}

eq.6.9 <- function(kappa,mu,zu,su){
    sigma <- exp(kappa)
    wu <- 1/(sigma^2+su^2)
    (1-sum((wu*(zu-mu))^2,na.rm=TRUE)/sum(wu,na.rm=TRUE))^2
}
