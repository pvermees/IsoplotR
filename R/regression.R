#' Linear regression of X,Y-variables with correlated errors
#'
#' Implements the unified regression algorithm of York et al. (2004)
#' which, although based on least squares, yields results that are
#' consistent with maximum likelihood estimates of Ludwig and
#' Titterington (1994)
#'
#' @param X vector of measurements
#' @param Y vector of measurements
#' @param sX standard errors of \code{X}
#' @param sY standard errors of \code{Y}
#' @param rXY correlation coefficients between X and Y
#' @return a five element list containing
#'
#' \code{a}: the intercept of the straight line fit
#'
#' \code{b}: the slope of the fit
#'
#' \code{sa}: the standard error of the intercept
#'
#' \code{sb}: the standard error of the slope
#' @references
#'
#' Ludwig, K. R., and D. M. Titterington. "Calculation of 230ThU
#' isochrons, ages, and errors." Geochimica et Cosmochimica Acta
#' 58.22 (1994): 5031-5042.
#' 
#' York, Derek, et al. "Unified equations for the slope,
#' intercept, and standard errors of the best straight line."
#' American Journal of Physics 72.3 (2004): 367-375.
#'
#' @importFrom stats lm
#' @examples
#'    X <- c(1.550,12.395,20.445,20.435,20.610,24.900,
#'           28.530,50.540,51.595,86.51,106.40,157.35)
#'    Y <- c(.7268,.7849,.8200,.8156,.8160,.8322,
#'           .8642,.9584,.9617,1.135,1.230,1.490)
#'    n <- length(X)
#'    sX <- X*0.01
#'    sY <- Y*0.005
#'    rXY <- rep(0.8,n)
#'    fit <- yorkfit(X,Y,sX,sY,rXY)
#'    covmat <- matrix(0,2,2)
#'    plot(range(X),fit$a+fit$b*range(X),type='l',ylim=range(Y))
#'    for (i in 1:n){
#'        covmat[1,1] <- sX[i]^2
#'        covmat[2,2] <- sY[i]^2
#'        covmat[1,2] <- rXY[i]*sX[i]*sY[i]
#'        covmat[2,1] <- covmat[1,2]
#'        ell <- ellipse(X[i],Y[i],covmat,alpha=0.05)
#'        polygon(ell)
#'    }
#' @export 
yorkfit <- function(X,Y,sX,sY,rXY){
    ab <- lm(Y ~ X)$coefficients # initial guess
    a <- ab[1]
    b <- ab[2]
    wX <- 1/sX^2
    wY <- 1/sY^2
    for (i in 1:50){ # 50 = maximum number of iterations
        bold <- b
        alpha <- sqrt(wX*wY)
        W <- wX*wY/(wX+b*b*wY-2*b*rXY*alpha)
        Xbar <- sum(W*X)/sum(W)
        Ybar <- sum(W*Y)/sum(W)
        U <- X-Xbar
        V <- Y-Ybar
        beta <- W*(U/wY+b*V/wX-(b*U+V)*rXY/alpha)
        b <- sum(W*beta*V)/sum(W*beta*U)
        if ((bold/b-1)^2 < 1e-15) break # convergence reached
    }
    a <- Ybar-b*Xbar
    x <- Xbar + beta
    xbar <- sum(W*x)/sum(W)
    u <- x-xbar
    sb <- sqrt(1/sum(W*u^2))
    sa <- sqrt(1/sum(W)+(xbar*sb)^2)
    out <- list()
    out$a <- a
    out$b <- b
    out$sa <- sa
    out$sb <- sb
    out
}

# get fitted x and y given a dataset X,Y,sX,sY,rXY,
# an intercept a and slope b. This function is useful
# for evaluating log-likelihoods of derived quantities
get.york.xy <- function(X,Y,sX,sY,rXY,a,b){
    wX <- 1/sX^2
    wY <- 1/sY^2
    alpha <- sqrt(wX*wY)
    W <- wX*wY/(wX+b*b*wY-2*b*rXY*alpha)
    Xbar <- sum(W*X)/sum(W)
    Ybar <- a + b*Xbar
    U <- X-Xbar
    V <- Y-Ybar
    beta <- W*(U/wY+b*V/wX-(b*U+V)*rXY/alpha)
    out <- cbind(Xbar+beta,Ybar+b*beta)
    out
}

UPb2york <- function(x,wetherill=TRUE){
    selection <- get.UPb.selection(wetherill)
    out <- list()
    nn <- nrow(x$x)
    out$X <- x$x[,selection[1]]
    out$Y <- x$x[,selection[2]]
    out$sX <- rep(0,nn)
    out$sY <- rep(0,nn)
    out$rXY <- rep(0,nn)
    for (i in 1:nn){
        covmat <- get.covmat.UPb(x,i)
        out$sX[i] <- sqrt(covmat[1,1])
        out$sY[i] <- sqrt(covmat[2,2])
        out$rXY[i] <- covmat[1,2]/(out$sX[i]*out$sY[i])
    }
    out
}
