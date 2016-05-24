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
#' @param dcu propagate the decay constant uncertainties?
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
concordia.age <- function(x,wetherill=TRUE,dcu=TRUE){
    X <- UPb.preprocess(x,wetherill)
    xy <- initialise.concordant.composition(x,wetherill)
    fit.comp <- optim(xy, LL.concordia.comp, x=X, method="BFGS", hessian=TRUE)
    out <- list()
    out$x <- fit.comp$par
    out$x.cov <- solve(fit.comp$hessian)
    selection <- names(X[[1]]$x)
    names(out$x) <- selection
    colnames(out$x.cov) <- selection
    rownames(out$x.cov) <- selection
    t.init <- initial.concordia.age.guess(out$x,wetherill)
    fit.age <- optim(t.init, LL.concordia.age, x=out$x, covmat=out$x.cov,
                     wetherill=wetherill, dcu=dcu, method="BFGS", hessian=TRUE)
    out$age <- fit.age$par
    out$age.err <- as.numeric(sqrt(solve(fit.age$hessian)))
    mswd <- mswd.concordia(X,out$x,out$x.cov,out$age,wetherill)
    out$mswd <- mswd$mswd
    out$p.value <- mswd$p.value
    out
}

initialise.concordant.composition <- function(x,wetherill=TRUE){
    selection <- get.UPb.selection(wetherill)
    colMeans(x$x[,selection])
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

mswd.concordia <- function(x,mu,covmat,age,wetherill=TRUE,dcu=TRUE){
    out <- list()
    SS.equivalence <- LL.concordia.comp(mu,x,TRUE)
    SS.concordance <- LL.concordia.age(age,mu,covmat,wetherill,mswd=TRUE,dcu=dcu)
    df.equivalence <- 2*length(x)-1
    df.concordance <- 1
    out$mswd <- list(equivalence = SS.equivalence/df.equivalence,
                     concordance = SS.concordance/df.concordance)
    out$p.value <- list(equivalence = 1-pchisq(SS.equivalence,df.equivalence),
                        concordance = 1-pchisq(SS.concordance,df.concordance))
    out
}

LL.concordia.comp <- function(mu,x,mswd=FALSE){
    out <- 0
    for (i in 1:length(x)){
        X <- matrix(x[[i]]$x-mu,1,2)
        covmat <- x[[i]]$cov
        if (mswd) out <- out + get.SS(X,covmat)
        else out <- out + LL.norm(X,covmat)
    }
    out
}

LL.concordia.age <- function(age,x,covmat,wetherill=TRUE,mswd=FALSE,dcu=TRUE){
    UPbratios <- get.ratios.UPb(age)
    selection <- get.UPb.labels(wetherill)
    X <- matrix(x[selection]-UPbratios$x[selection],1,2)
    COVMAT <- covmat[selection,selection]
    if (dcu) COVMAT <- COVMAT + UPbratios$cov[selection,selection]
    if (mswd) out <- get.SS(X,COVMAT)
    else out <- LL.norm(X,COVMAT)
    out
}

get.UPb.labels <- function(wetherill=TRUE){
    if (wetherill) selection <- c('Pb207U235','Pb206U238')
    else selection <- c('U238Pb206','Pb207Pb206')
    selection
}

#' Linear regression on a U-Pb concordia diagram
#'
#' Performs linear regression of U-Pb data on Wetherill and
#' Tera-Wasserburg concordia diagrams. Computes the upper and lower
#' intercept ages (for Wetherill) or the lower intercept age and the
#' \eqn{^{207}}Pb/\eqn{^{206}}Pb intercept (for Tera-Wasserburg),
#' taking into account error correlations and decay constant
#' uncertainties.
#'
#' @param x an object of class \code{UPb}
#' @param wetherill boolean flag to indicate whether the data should
#'     be evaluated in Wetherill (\code{TRUE}) or Tera-Wasserburg
#'     (\code{FALSE}) space
#' @param dcu propagate the decay constant uncertainties?
#' @return a list with the following items:
#'
#' \code{x}: a two element vector with the upper and lower intercept
#' ages (if wetherill==TRUE) or the lower intercept age and
#' \eqn{^{207}}Pb/\eqn{^{206}}Pb intercept (for Tera-Wasserburg)
#'
#' \code{cov}: the covariance matrix of the elements in \code{x}
#'
#' @importFrom stats optim optimHess
#' @examples
#' data(UPb)
#' fit <- discordia.age(UPb)
#' print(paste('lower intercept = ',fit$x[1],'+/-',sqrt(fit$cov[1,1]),'Ma'))
#' @export
discordia.age <- function(x,wetherill=TRUE,dcu=TRUE){
    out <- list()
    d <- UPb2york(x,wetherill)
    X <- UPb.preprocess(x,wetherill)
    fit <- yorkfit(d$X,d$Y,d$sX,d$sY,d$rXY)
    itt <- concordia.intersection(fit,wetherill)
    hess <- optimHess(itt,LL.concordia.intersection,
                      d=d,X=X,wetherill=wetherill,dcu=dcu)
    out$x <- itt
    out$cov <- solve(hess)
    out
}

# itt = output of the yorkfit function
# d = output of UPb2york
# X = output of UPb.preprocess
LL.concordia.intersection <- function(itt,d,X,wetherill,dcu){
    LL <- 0
    selection <- get.UPb.labels(wetherill)
    XYl <- get.ratios.UPb(itt[1])$x[selection]
    if (wetherill){
        XYu <- get.ratios.UPb(itt[2])$x[selection]
        b <- (XYu[2]-XYl[2])/(XYu[1]-XYl[1])
        a <- XYu[2]-b*XYu[1]
    } else {
        a <- itt[2]
        b <- (XYl[2]-a)/XYl[1]
    }
    xy <- get.york.xy(d$X,d$Y,d$sX,d$sY,d$rXY,a,b)
    mu <- cbind(d$X,d$Y)-xy
    if (wetherill) disc <- (XYu[1]-xy[,1])/(XYu[1]-XYl[1])
    else disc <- (XYl[1]-xy[,1])/XYl[1]
    for (i in 1:length(X)){
        COVMAT <- X[[i]]$cov
        if (dcu) {
            dcomp <- discordant.composition(disc[i],itt[1],itt[2],wetherill)
            COVMAT <- COVMAT + dcomp$cov
        }
        LL <- LL + LL.norm(matrix(mu[i,],1,2),COVMAT)        
    }
    LL
}

# returns the lower and upper intercept age (for Wetherill concordia)
# or the lower intercept age and 207Pb/206Pb intercept (for Tera-Wasserburg)
concordia.intersection <- function(fit,wetherill){
    if (wetherill){
        search.range <- c(0,10000)
        midpoint <- stats::optimize(intersection.misfit, search.range,
                                    a=fit$a, b=fit$b, wetherill=wetherill)$minimum
        range1 <- c(-1000,midpoint)
        range2 <- c(midpoint,10000)
        tt <- search.range # tl, tu
        tt[1] <- stats::uniroot(intersection.misfit, range1, 
                                a=fit$a, b=fit$b, wetherill=wetherill)$root
        tt[2] <- stats::uniroot(intersection.misfit, range2, 
                                a=fit$a, b=fit$b, wetherill=wetherill)$root
        out <- tt
    } else {
        search.range <- c(1/10000,10000)
        it <- c(1,fit$a) # tl, 7/6 intercept
        if (fit$b<0) { # negative slope => two intersections with concordia line
            midpoint <- stats::optimize(intersection.misfit, search.range,
                                        a=fit$a, b=fit$b, wetherill=wetherill)$minimum
            search.range[2] <- midpoint
            it[1] <- stats::uniroot(intersection.misfit, search.range, 
                                    a=fit$a, b=fit$b, wetherill=wetherill)$root
            out <- it
        } else {
            it[1] <- stats::uniroot(intersection.misfit, search.range,
                                    a=fit$a, b=fit$b, wetherill=wetherill)$root
        }
        out <- it
    }
    out
}

# returns misfit of a proposed age and the intersection between the
# discordia and concordia lines
intersection.misfit <- function(age,a,b,wetherill){
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    R <- I.R('U238U235')[1]
    if (wetherill){
        out <- a-b+1 + b*exp(l5*age) - exp(l8*age)
    } else {
        out <- (exp(l5*age)-1)/(exp(l8*age)-1) - a*R - b*R/(exp(l8*age)-1)
    }
    out
}

# find the composition of a point that is 100d% discordant
# between upper intercept and lower intercept age
# (if wetherill=TRUE) or between the 207/206 intercept
# and the lower intercept age (if wetherill=FALSE)
discordant.composition <- function(d,tl,itu,wetherill=TRUE){
    out <- list()
    l5 <- lambda('U235')[1]
    l8 <- lambda('U238')[1]
    R <- I.R('U238U235')[1]
    if (wetherill){
        X <- d*(exp(l5*tl)-1) + (1-d)*(exp(l5*itu)-1)
        Y <- d*(exp(l8*tl)-1) + (1-d)*(exp(l8*itu)-1)
        dXdl5 <- d*tl*exp(l5*tl) + (1-d)*itu*exp(l5*itu)
        dXdl8 <- 0
        dXdR <- 0
        dYdl5 <- 0
        dYdl8 <- d*tl*exp(l8*tl) + (1-d)*itu*exp(l8*itu)
        dYdR <- 0
    } else {
        X <- (1-d)/(exp(l8*tl)-1)
        Y <- d*itu + ((1-d)/R)*(exp(l5*tl)-1)/(exp(l8*tl)-1)
        dXdl5 <- 0
        dXdl8 <- -(1-d)*tl*exp(l8*tl)/(exp(l8*tl)-1)^2
        dXdR <- 0
        dYdl5 <- ((1-d)/R)*tl*exp(l5*tl)/(exp(l8*tl)-1)
        dYdl8 <- -((1-d)/R)*tl*exp(l8*tl)*(exp(l5*tl)-1)/(exp(l8*tl)-1)^2
        dYdR <- -((1-d)/R^2)*(exp(l5*tl)-1)/(exp(l8*tl)-1)
    }
    J <- matrix(0,2,3)
    J[1,1] <- dXdl5
    J[1,2] <- dXdl8
    J[1,3] <- dXdR
    J[2,1] <- dYdl5
    J[2,2] <- dYdl8
    J[2,3] <- dYdR
    E <- matrix(0,3,3)
    E[1,1] <- lambda('U235')[2]^2
    E[2,2] <- lambda('U238')[2]^2
    E[3,3] <- I.R('U238U235')[2]^2
    covmat <- J %*% E %*% t(J)
    out$x <- c(X,Y)
    out$cov <- covmat
    labels <- get.UPb.labels(wetherill)
    names(out$x) <- labels
    rownames(out$cov) <- labels
    colnames(out$cov) <- labels
    out
}

# negative multivariate log likelihood to be fed into R's optim function
LL.norm <- function(x,covmat){
    log(2*pi) + 0.5*determinant(covmat,logarithmic=TRUE)$modulus + 0.5*get.SS(x,covmat)
}

get.SS <- function(x,covmat){
    x %*% solve(covmat) %*% t(x)
}

discordia.plot <- function(fit,wetherill=TRUE){
    selection <- get.UPb.labels(wetherill)
    comp.l <- get.ratios.UPb(fit$x[1])$x[selection]
    if (wetherill) {
        comp.u <- get.ratios.UPb(fit$x[2])$x[selection]
    } else {
        comp.u <- c(0,fit$x[2])
    }
    X <- c(comp.l[1],comp.u[1])
    Y <- c(comp.l[2],comp.u[2])
    graphics::lines(X,Y)
}
