#' Calculate U-Th-He (and fission track) central ages and compositions
#'
#' Computes the geometric mean composition of a set of fission track
#' or U-Th-He data and returns the corresponding age and fitting
#' parameters.
#'
#' @param x an object of class \code{UThHe} or \code{fissiontracks}
#' @param ... optional arguments
#' @return a list containing the following items:
#'
#' \describe{
#'
#' \item{mswd}{the reduced Chi-square statistic of data concordance,
#' i.e. \eqn{mswd=SS/(2n-2)}, where \eqn{SS} is the sum of squares
#' of the log[U/He]-log[Th/He] compositions and \eqn{n} is the number
#' of samples. If \code{x} has class \code{fissiontracks}, then
#' \eqn{mswd=X^2/(n-1)}, where \eqn{X^2} is a Chi-square statistic of
#' the EDM data or ICP ages.}
#'
#' \item{p.value}{the p-value of a Chi-square test with \eqn{n-2}
#' degrees of freedom}
#'
#' \item{age}{a two-column vector with the central age and its
#' standard error.}
#'
#' }
#'
#' Additionally, if \code{x} has class \code{UThHe}:
#'
#' \describe{
#'
#' \item{uvw}{(if the input data table contains Sm) or \strong{uv} (if
#' it doesn't): the geometric mean log[U/He], log[Th/He] (,
#' log[Sm/He]) and log[Sm/He] composition}
#'
#' \item{covmat}{the covariance matrix of \code{uvw} or \code{uv}}
#'
#' }
#'
#' OR, if \code{x} has class \code{fissiontracks}:
#'
#' \describe{
#'
#' \item{disp}{the (over)dispersion of the ages (value between 0 and
#' 1)}
#'
#' }
#'
#' @examples
#' data(examples)
#' print(central(examples$UThHe)$age)
#' @rdname central
#' @export
central <- function(x,...){ UseMethod("central",x) }
#' @rdname central
#' @export
central.default <- function(x,...){
    stop('Invalid input into central(...) function');
}
#' @rdname central
#' @export
central.UThHe <- function(x,...){
    out <- list()
    ns <- nrow(x)
    df <- 2*(ns-1)
    doSm <- doSm(x)
    if (doSm){
        uvw <- UThHe2uvw(x)
        fit <- stats::optim(c(0,0,0),SS.UThHe.uvw,method='BFGS',
                            hessian=TRUE,x=x)
        out$uvw <- fit$par
        out$covmat <- solve(fit$hessian)
        nms <- c('u','v','w')
        cc <- uvw2UThHe(out$uvw,out$covmat)
        out$age <- get.UThHe.age(cc['U'],cc['sU'],cc['Th'],cc['sTh'],
                                 cc['He'],cc['sHe'],cc['Sm'],cc['sSm'])
        SS <- SS.UThHe.uv(out$uvw[1:2],x)
    } else {
        uv <- UThHe2uv(x)
        fit <- stats::optim(c(0,0),SS.UThHe.uv,method='BFGS',
                            hessian=TRUE,x=x)
        out$uv <- fit$par
        out$covmat <- solve(fit$hessian)
        nms <- c('u','v')
        cc <- uv2UThHe(out$uv,out$covmat)
        out$age <- get.UThHe.age(cc['U'],cc['sU'],
                                 cc['Th'],cc['sTh'],
                                 cc['He'],cc['sHe'])
        SS <- SS.UThHe.uv(out$uv,x)
    }
    names(out$uv) <- nms
    colnames(out$covmat) <- nms
    rownames(out$covmat) <- nms
    out$mswd <- SS/df
    out$p.value <- 1-pchisq(SS,df)
    out
}
#' @param mineral setting this parameter to either \code{apatite} or
#'     \code{zircon} changes the default efficiency factor, initial
#'     fission track length and density to preset values (only affects
#'     results if \code{x$format=2}.)
#' @rdname central
#' @export
central.fissiontracks <- function(x,mineral=NA,...){
    out <- list()
    L8 <- lambda('U238')[1]
    sigma <- 0.15 # convenient starting value
    if (x$format<2){
        Nsj <- x$x[,'Ns']
        Ns <- sum(Nsj)
        Nij <- x$x[,'Ni']
        Ni <- sum(Nij)
        num <- (Nsj*Ni-Nij*Ns)^2
        den <- Nsj+Nij
        df <- length(Nsj)-1
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
    } else if (x$format>1){
        tst <- age(x,exterr=FALSE,mineral=mineral)
        zu <- log(tst[,'t'])
        su <- tst[,'s[t]']/tst[,'t']
        for (i in 1:30){ # page 100 of Galbraith (2005)
            wu <- 1/(sigma^2+su^2)
            mu <- sum(wu*zu,na.rm=TRUE)/sum(wu,na.rm=TRUE)
            fit <- stats::optimize(eq.6.9,c(0,10),mu=mu,zu=zu,su=su)
            sigma <- fit$minimum
        }
        tt <- exp(mu)
        st <- 1/sqrt(sum(wu,na.rm=TRUE))
        if (x$format==2)
            st <- tt * sqrt((st/tt)^2 + (x$zeta[2]/x$zeta[1])^2)
        df <- length(zu)-1
        Chi2 <- sum((zu/su)^2,na.rm=TRUE)-(sum(zu/su^2,na.rm=TRUE)^2)/sum(1/su^2,na.rm=TRUE)
    }
    out$age <- c(tt,st)
    out$disp <- sigma
    out$mswd <- Chi2/df
    out$p.value <- 1-pchisq(Chi2,df)
    out
}

eq.6.9 <- function(sigma,mu,zu,su){
    wu <- 1/(sigma^2+su^2)
    (1-sum((wu*(zu-mu))^2,na.rm=TRUE)/sum(wu,na.rm=TRUE))^2
}
