#' Linear regression of X,Y-variables with correlated errors
#'
#' Implements the unified regression algorithm of York et al. (2004)
#' which, although based on least squares, yields results that are
#' consistent with maximum likelihood estimates of Ludwig and
#' Titterington (1994)
#'
#' @param X EITHER a vector with the X-values OR a 5-column matrix
#'     with the X-values, the analytical uncertainties of the
#'     X-values, the Y-values, the analytical uncertainties of the
#'     Y-values, and the correlation coefficients of the X- and
#'     Y-values.
#' @param sX standard errors of \code{X} OR \code{NULL} if \code{X} is a
#'     matrix or data frame
#' @param Y vector of measurements with the same length as \code{X} OR
#'     \code{NULL} if \code{X} is a matrix or data frame
#' @param sY standard errors of \code{Y} OR \code{NULL} if \code{X} is a
#'     matrix or data frame
#' @param rXY correlation coefficients between \code{X} and \code{Y}
#'     OR \code{NULL} if \code{X} is a matrix or data frame
#' @return a three-element list of vectors containing:
#'     \describe{
#'     \item{a}{the intercept of the straight line fit and its
#'     standard error}
#'     \item{b}{the slope of the fit and its standard error}
#'     \item{cov.ab}{the covariance of the slope and intercept}
#'     }
#' @references
#'
#' Ludwig, K. R., and D. M. Titterington. "Calculation of \eqn{^{230}}Th/U
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
#'    fit <- yorkfit(X,sX,Y,sY,rXY)
#'    covmat <- matrix(0,2,2)
#'    plot(range(X),fit$a[1]+fit$b[1]*range(X),type='l',ylim=range(Y))
#'    for (i in 1:n){
#'        covmat[1,1] <- sX[i]^2
#'        covmat[2,2] <- sY[i]^2
#'        covmat[1,2] <- rXY[i]*sX[i]*sY[i]
#'        covmat[2,1] <- covmat[1,2]
#'        ell <- ellipse(X[i],Y[i],covmat,alpha=0.05)
#'        polygon(ell)
#'    }
#' @export
yorkfit <- function(X,sX=NULL,Y=NULL,sY=NULL,rXY=NULL){
    if (hasClass(X,'matrix') | hasClass(X,'data.frame'))
        return(yorkfit(X[,1],X[,2],X[,3],X[,4],X[,5]))
    ab <- lm(Y ~ X)$coefficients # initial guess
    a <- ab[1]
    b <- ab[2]
    wX <- 1/sX^2
    wY <- 1/sY^2
    for (i in 1:50){ # 50 = maximum number of iterations
        bold <- b
        alpha <- sqrt(wX*wY)
        W <- wX*wY/(wX+b*b*wY-2*b*rXY*alpha)
        Xbar <- sum(W*X,na.rm=TRUE)/sum(W,na.rm=TRUE)
        Ybar <- sum(W*Y,na.rm=TRUE)/sum(W,na.rm=TRUE)
        U <- X-Xbar
        V <- Y-Ybar
        beta <- W*(U/wY+b*V/wX-(b*U+V)*rXY/alpha)
        b <- sum(W*beta*V,na.rm=TRUE)/sum(W*beta*U,na.rm=TRUE)
        if ((bold/b-1)^2 < 1e-15) break # convergence reached
    }
    a <- Ybar-b*Xbar
    x <- Xbar + beta
    xbar <- sum(W*x,na.rm=TRUE)/sum(W,na.rm=TRUE)
    u <- x-xbar
    sb <- sqrt(1/sum(W*u^2,na.rm=TRUE))
    sa <- sqrt(1/sum(W,na.rm=TRUE)+(xbar*sb)^2)
    out <- list()
    out$a <- c(a,sa)
    out$b <- c(b,sb)
    out$cov.ab <- -Xbar*sb^2
    out$mswd <- 0
    mswd <- get.york.mswd(X,sX,Y,sY,rXY,a,b)
   out$mswd <- mswd$mswd
    out$p.value <- mswd$p.value
    out
}

get.york.mswd <- function(X,sX,Y,sY,rXY,a,b){
    xy <- get.york.xy(X,sX,Y,sY,rXY,a,b)
    X2 <- 0
    nn <- length(X)
    for (i in 1:nn){
        E <- cor2cov(sX[i],sY[i],rXY[i])
        x <- matrix(c(X[i]-xy[i,1],Y[i]-xy[i,2]),1,2)
        if (!any(is.na(x)))
            X2 <- X2 + 0.5*x %*% solve(E) %*% t(x)
    }
    df <- (2*nn-2)
    out <- list()
    out$mswd <- X2/df
    out$p.value <- 1-pchisq(X2,df)
    out
}

# get fitted x and y given a dataset X,sX,Y,sY,rXY,
# an intercept a and slope b. This function is useful
# for evaluating log-likelihoods of derived quantities
get.york.xy <- function(X,sX,Y,sY,rXY,a,b){
    wX <- 1/sX^2
    wY <- 1/sY^2
    alpha <- sqrt(wX*wY)
    W <- wX*wY/(wX+b*b*wY-2*b*rXY*alpha)
    Xbar <- sum(W*X,na.rm=TRUE)/sum(W,na.rm=TRUE)
    Ybar <- a + b*Xbar
    U <- X-Xbar
    V <- Y-Ybar
    beta <- W*(U/wY+b*V/wX-(b*U+V)*rXY/alpha)
    out <- cbind(Xbar+beta,Ybar+b*beta)
    out
}

data2york <- function(x,selection=NA){
    out <- list()
    nn <- nrow(x$x)
    if (any(is.na(selection)))
        selection <- c(1,2)
    out$X <- x$x[,selection[1]]
    out$Y <- x$x[,selection[2]]
    out$sX <- rep(0,nn)
    out$sY <- rep(0,nn)
    out$rXY <- rep(0,nn)
    for (i in 1:nn){
        covmat <- get.covmat(x,i)[selection,selection]
        out$sX[i] <- sqrt(covmat[1,1])
        out$sY[i] <- sqrt(covmat[2,2])
        out$rXY[i] <- stats::cov2cor(covmat)[1,2]
    }
    out
}
