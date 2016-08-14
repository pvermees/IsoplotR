#' Calculate U-Th-He (and fission track) central ages and compositions
#'
#' Computes the geometric mean composition of a set of fission track
#' or U-Th-He data and returns the corresponding age and fitting
#' parameters.
#'
#' @param x an object of class \code{UThHe} or \code{fissiontracks}
#' @param ... optional arguments
#' @return if \code{x} has class \code{UThHe}, retuns a list
#'     containing the following items:
#'
#' \code{uvw} (if the input data table contains Sm) or \code{uv} (if it doesn't):
#' the geometric mean log[U/He], log[Th/He] (, log[Sm/He]) and log[Sm/He]
#' composition
#'
#' \code{covmat}: the covariance matrix of \code{uvw} or \code{uv}
#'
#' \code{mswd}: the reduced Chi-square statistic of data concordance,
#' i.e. mswd = SS/(2n-2) where SS is the sum of squares of the
#' log[U/He]-log[Th/He] compositions and n is the number of samples.
#'
#' \code{p.value}: the p-value of a Chi-square test with n-2 degrees
#' of freedom
#'
#' \code{age}: a two-column vector with the central age and its
#' standard error.
#'
#' @examples
#' data(examples)
#' print(central(examples$UThHe)$age)
#' @rdname central
#' @export
central <- function(x,...){ UseMethod("central",x) }
#' @rdname central
#' @export
central.default <- function(x,...){ stop('Invalid input into central(...) function'); }
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
