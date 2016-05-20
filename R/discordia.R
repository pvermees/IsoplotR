#' Calculate U-Pb concordia ages
#'
#' Evaluates the equivalence of multiple
#' (\eqn{^{206}}Pb/\eqn{^{238}}U-\eqn{^{207}}Pb/\eqn{^{235}}U or
#' \eqn{^{207}}Pb/\eqn{^{206}}Pb-\eqn{^{206}}Pb/\eqn{^{238}}U)
#' compositions, computes the weighted mean isotopic composition and
#' the corresponding concordia age using the method of maximum
#' likelihood, computes the mswd of equivalence and concordance and
#' their respective Chi-squared p-values.
#'
#' @param x an object of class \code{UPb}
#' @param wetherill boolean flag to indicate whether the data should
#'     be evaluated in Wetherill (\code{TRUE}) or Tera-Wasserburg
#'     (\code{FALSE}) space
#' @return a list with the following items:
#'
#' \code{x}: a named vector with the weighted mean U-Pb composition
#'
#' \code{x.cov}: the covariance matrix of the mean U-Pb composition
#'
#' \code{age}: the concordia age (in Ma)
#'
#' \code{age.err}: the standard error of the concordia age
#'
#' \code{mswd}: a list with two items (\code{equivalence} and
#' \code{concordance}) containing the MSWD (Mean of the Squared
#' Weighted Deviates, a.k.a the reduced Chi-squared statistic outside
#' of geochronology) of isotopic equivalence and age concordance,
#' respectively.
#'
#' \code{p.value}: a list with two items (\code{equivalence} and
#' \code{concordance}) containing the p-value of the Chi-square test
#' for isotopic equivalence and age concordance, respectively.
#' @importFrom stats optim
#' @examples
#' data(UPb)
#' fit <- concordia.age(UPb)
#' print(paste('age = ',fit$age,'+/-',fit$age.err,'Ma, MSWD = ',fit$mswd))
#' @export
concordia.age <- function(x,wetherill=TRUE){
    X <- UPb.preprocess(x,wetherill)
    fit.comp <- optim(c(1,1), LL.concordia, x=X, method="BFGS", hessian=TRUE)
    out <- list()
    out$x <- fit.comp$par
    out$x.cov <- solve(fit.comp$hessian)
    selection <- names(X[[1]]$x)
    names(out$x) <- selection
    colnames(out$x.cov) <- selection
    rownames(out$x.cov) <- selection
    t.init <- initial.concordia.age.guess(out$x,wetherill)
    fit.age <- optim(t.init, LL.concordia.age, x=out$x, covmat=out$x.cov,
                     wetherill=wetherill, method="BFGS", hessian=TRUE)
    out$age <- fit.age$par
    out$age.err <- as.numeric(sqrt(solve(fit.age$hessian)))
    mswd <- mswd.concordia(X,out$x,out$x.cov,out$age,wetherill)
    out$mswd <- mswd$mswd
    out$p.value <- mswd$p.value
    out
}

initial.concordia.age.guess <- function(x,wetherill=TRUE){
    if (wetherill){
        Pb207U235age <- get.Pb207U235age(x['Pb207U235'])
        Pb206U238age <- get.Pb206U238age(x['Pb206U238'])
        out <- mean(c(Pb207U235age,Pb206U238age))
    } else {
        Pb206U238age <- get.Pb206U238age(1/x['U238Pb206'])
        Pb207Pb206age <- get.Pb207Pb206age(x['Pb207Pb206'])
        out <- mean(c(Pb206U238age,Pb207Pb206age))
    }
    out
}

get.UPb.selection <- function(wetherill=TRUE){
    if (wetherill) selection <- c('Pb207U235','Pb206U238')
    else selection <- c('U238Pb206','Pb207Pb206')
    selection
}

UPb.preprocess <- function(x,wetherill){
    selection <- get.UPb.selection(wetherill)
    out <- list()
    for (i in 1:nrow(x$x)){
        X <- x$x[i,selection]
        covmat <- get.covmat.UPb(x,i)[selection,selection]
        out[[i]] <- list(x=X, cov=covmat)
    }
    out   
}

mswd.concordia <- function(x,mu,covmat,age,wetherill=TRUE){
    out <- list()
    SS.equivalence <- LL.concordia(mu,x,TRUE)
    SS.concordance <- LL.concordia.age(age,mu,covmat,wetherill,TRUE)
    df.equivalence <- 2*length(x)-1
    df.concordance <- 1
    out$mswd <- list(equivalence = SS.equivalence/df.equivalence,
                     concordance = SS.concordance/df.concordance)
    out$p.value <- list(equivalence = 1-pchisq(SS.equivalence,df.equivalence),
                        concordance = 1-pchisq(SS.concordance,df.concordance))
    out
}

LL.concordia <- function(mu,x,mswd=FALSE){
    out <- 0
    for (i in 1:length(x)){
        X <- matrix(x[[i]]$x-mu,1,2)
        covmat <- x[[i]]$cov
        if (mswd) out <- out + get.SS(X,covmat)
        else out <- out - LL.norm(X,covmat)
    }
    out
}

LL.concordia.age <- function(age,x,covmat,wetherill=TRUE,mswd=FALSE){
    UPbratios <- get.ratios.UPb(age)
    if (wetherill) selection <- c('Pb207U235','Pb206U238')
    else selection <- c('U238Pb206','Pb207Pb206')
    X <- matrix(x[selection]-UPbratios$x[selection],1,2)
    COVMAT <- UPbratios$cov[selection,selection] + covmat[selection,selection]
    if (mswd) out <- get.SS(X,COVMAT)
    else out <- -LL.norm(X,COVMAT)
    out
}

LL.norm <- function(x,covmat){
    - log(2*pi) - 0.5*determinant(covmat,logarithmic=TRUE)$modulus - 0.5*get.SS(x,covmat)
                  
}

get.SS <- function(x,covmat){
    x %*% solve(covmat) %*% t(x)
}
