#' Calculate U-Th-He and fission track central ages and compositions
#'
#' Computes the geometric mean composition of a set of fission track
#' or U-Th-He data and returns the corresponding age and fitting
#' parameters.
#'
#' @param x an object of class \code{UThHe} or \code{fissiontracks},
#'     OR a 2-column matrix with (strictly positive) values and
#'     uncertainties
#' @param alpha cutoff value for confidence intervals
#' @param ... optional arguments
#' @return if \code{x} has class \code{UThHe}, a list containing the
#'     following items:
#'
#' \describe{
#' \item{age}{a four-element vector with:
#'
#' \code{t}: the central age
#'
#' \code{s[t]}: the standard error of \code{s[t]}
#'
#' \code{ci[t]}: the \eqn{100(1-\alpha/2)\%} confidence interval for
#' \code{t} for the appropriate number of degrees of freedom
#'
#' \code{disp[t]}: the \eqn{100(1-\alpha/2)\%} confidence interval
#' enhanced by a factor of \eqn{\sqrt{MSWD}}.
#' }
#'
#' \item{mswd}{the reduced Chi-square statistic of data concordance,
#' i.e. \eqn{mswd=SS/df}, where \eqn{SS} is the sum of squares of
#' the log[U/He]-log[Th/He] compositions}
#'
#' \item{df}{the degrees of freedom (\eqn{2n-2})}
#'
#' \item{p.value}{the p-value of a Chi-square test with \code{df}
#' degrees of freedom}
#'
#' \item{uvw}{(if the input data table contains Sm) or \strong{uv} (if
#' it doesn't): the geometric mean log[U/He], log[Th/He] (,
#' log[Sm/He]) and log[Sm/He] composition}
#'
#' \item{covmat}{the covariance matrix of \code{uvw} or \code{uv}}
#'
#' \item{tfact}{the \eqn{100(1-\alpha/2)\%} percentile of the t-
#' distribution for \code{df} degrees of freedom}
#'
#' }
#'
#' OR, otherwise:
#'
#' \describe{
#'
#' \item{age}{a three-element vector with:
#'
#' \code{t}: the central age
#'
#' \code{s[t]}: the standard error of \code{s[t]}
#'
#' \code{ci[t]}: the \eqn{100(1-\alpha/2)\%} confidence interval for
#' \code{t} for the appropriate number of degrees of freedom }
#'
#' \item{disp}{a two-element vector with the overdispersion (standard
#' deviation) of the excess scatter, and the corresponding
#' \eqn{100(1-\alpha/2)\%} confidence interval for the appropriate
#' degrees of freedom.}
#'
#' \item{mswd}{the reduced Chi-square statistic of data concordance,
#' i.e. \eqn{mswd=X^2/df}, where \eqn{X^2} is a Chi-square statistic
#' of the EDM data or ages}
#'
#' \item{df}{the degrees of freedom (\eqn{n-2})}
#'
#' \item{p.value}{the p-value of a Chi-square test with \code{df}
#' degrees of freedom}
#'
#' }
#'
#' @examples
#' data(examples)
#' print(central(examples$UThHe)$age)
#' 
#' @references Galbraith, R.F. and Laslett, G.M., 1993. Statistical
#'     models for mixed fission track ages. Nuclear tracks and
#'     radiation measurements, 21(4), pp.459-470.
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
    sigma <- 0.15 # convenient starting value
    zu <- log(x[,1])
    su <- x[,2]/x[,1]
    for (i in 1:30){ # page 100 of Galbraith (2005)
        wu <- 1/(sigma^2+su^2)
        mu <- sum(wu*zu,na.rm=TRUE)/sum(wu,na.rm=TRUE)
        fit <- stats::optimize(eq.6.9,c(0,10),mu=mu,zu=zu,su=su)
        sigma <- fit$minimum
    }
    tt <- exp(mu)
    st <- tt/sqrt(sum(wu,na.rm=TRUE))
    Chi2 <- sum((zu/su)^2,na.rm=TRUE)-(sum(zu/su^2,na.rm=TRUE)^2)/
        sum(1/su^2,na.rm=TRUE)
    out <- list()
    # remove two d.o.f. for mu and sigma
    out$df <- length(zu)-2
    # add back one d.o.f. for the homogeneity test
    out$mswd <- Chi2/(out$df+1)
    out$p.value <- 1-stats::pchisq(Chi2,out$df+1)
    out$age <- c(tt,st,stats::qt(1-alpha/2,out$df)*st)
    out$disp <- c(sigma,stats::qnorm(1-alpha/2)*sigma)
    names(out$age) <- c('t','s[t]','ci[t]')
    names(out$disp) <- c('s','ci')
    out
}
#' @rdname central
#' @export
central.UThHe <- function(x,alpha=0.05,...){
    out <- list()
    ns <- nrow(x)
    doSm <- doSm(x)
    out$age <- rep(NA,4)
    names(out$age) <- c('t','s[t]','ci[t]','disp[t]')
    if (doSm){
        uvw <- UThHe2uvw(x)
        fit <- stats::optim(c(0,0,0),SS.UThHe.uvw,method='BFGS',
                            hessian=TRUE,x=x)
        out$uvw <- fit$par
        out$covmat <- solve(fit$hessian)
        nms <- c('u','v','w')
        cc <- uvw2UThHe(out$uvw,out$covmat)
        out$age[c('t','s[t]')] <- get.UThHe.age(cc['U'],cc['sU'],
                                                cc['Th'],cc['sTh'],
                                                cc['He'],cc['sHe'],
                                                cc['Sm'],cc['sSm'])
    } else {
        uv <- UThHe2uv(x)
        fit <- stats::optim(c(0,0),SS.UThHe.uv,method='BFGS',
                            hessian=TRUE,x=x)
        out$uvw <- fit$par
        SS <- SS.UThHe.uv(out$uv[1:2],x)
        out$covmat <- solve(fit$hessian)
        nms <- c('u','v')
        cc <- uv2UThHe(out$uvw,out$covmat)
        out$age[c('t','s[t]')] <-
            get.UThHe.age(cc['U'],cc['sU'],
                          cc['Th'],cc['sTh'],
                          cc['He'],cc['sHe'])
    }
    # the overdispersion calculation does not use Sm
    SS <- SS.UThHe.uv(out$uvw[1:2],x)
    out$df <- 2*(ns-1)
    out$mswd <- SS/out$df
    out$tfact <- stats::qt(1-alpha/2,out$df)
    if (doSm){
        cco <- uvw2UThHe(out$uvw,out$mswd*out$covmat)
        out$age['disp[t]'] <-
            out$tfact*get.UThHe.age(cco['U'],cco['sU'],
                                    cco['Th'],cco['sTh'],
                                    cco['He'],cco['sHe'],
                                    cco['Sm'],cco['sSm'])[2]
    } else {
        cco <- uv2UThHe(out$uvw,out$mswd*out$covmat)
        out$age['disp[t]'] <-
            out$tfact*get.UThHe.age(cco['U'],cco['sU'],
                                    cco['Th'],cco['sTh'],
                                    cco['He'],cco['sHe'])[2]
    }
    out$age['ci[t]'] <- out$tfact*out$age['s[t]']
    names(out$uvw) <- nms
    colnames(out$covmat) <- nms
    rownames(out$covmat) <- nms
    out$p.value <- 1-stats::pchisq(SS,out$df)
    out
}
#' @param mineral setting this parameter to either \code{apatite} or
#'     \code{zircon} changes the default efficiency factor, initial
#'     fission track length and density to preset values (only affects
#'     results if \code{x$format=2}.)
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
        # remove two d.o.f. for mu and sigma
        out$df <- length(Nsj)-2
        # add back one d.o.f. for homogeneity test
        out$mswd <- Chi2/(out$df+1)
        out$p.value <- 1-stats::pchisq(Chi2,out$df+1)
        out$age <- c(tt,st,stats::qt(1-alpha/2,out$df)*st)        
        out$disp <- c(sigma,stats::qnorm(1-alpha/2)*sigma)
        names(out$age) <- c('t','s[t]','ci[t]')
        names(out$disp) <- c('s','ci')
    } else if (x$format>1){
        tst <- age(x,exterr=FALSE,mineral=mineral)
        out <- central.default(tst,alpha=alpha)
    }
    out
}

eq.6.9 <- function(sigma,mu,zu,su){
    wu <- 1/(sigma^2+su^2)
    (1-sum((wu*(zu-mu))^2,na.rm=TRUE)/sum(wu,na.rm=TRUE))^2
}
