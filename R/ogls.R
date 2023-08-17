#' @title
#' Omnivariant Generalised Least-Squares Regression
#'
#' @description
#' Linear regression with inter-sample error correlations.
#'
#' @param x either a \code{n x (n+1)} matrix obtained by prepending a
#'     vector of alternating \code{X,Y}-values to its covariance
#'     matrix OR an object of class \code{other} with
#'     \code{x$format=6}.
#' 
#' @param ... optional arguments
#' @return a list of the slope and intercept of the best fit line as
#'     well as their standard errors and covariance.
#' @examples
#' fn <- system.file('UW137.csv',package='IsoplotR')
#' UW137 <- read.data(fn,method='other',format=6)
#' fit <- ogls(UW137)
#' @rdname ogls
#' @export
ogls <- function(x,...){ UseMethod("ogls",x) }
#' @rdname ogls
#' @export
ogls.default <- function(x,random.effects=NA,...){
    out <- list()
    ydat <- data2york(x=x,format=6)
    yfit <- york(ydat)
    init <- c(yfit$a[1],yfit$b[1])
    omega <- solve(x[,-1])
    fit <- stats::optim(init,LL.ogls,dat=x,omega=omega,hessian=TRUE)
    covmat <- solve(fit$hessian)
    out$a <- c('a'=unname(fit$par[1]),'s[a]'=unname(sqrt(covmat[1,1])))
    out$b <- c('b'=unname(fit$par[2]),'s[b]'=unname(sqrt(covmat[2,2])))
    out$cov.ab <- covmat[1,2]
    out$type <- "ogls"
    out <- append(out,getMSWD(X2=fit$value,df=nrow(x)/2-2))
    out
}
#' @rdname ogls
#' @export
ogls.other <- function(x,random.effects=NA,...){
    ogls(x=x$x,random.effects=random.effects)
}

LL.ogls <- function(ab,dat,omega){
    a <- ab[1]
    b <- ab[2]
    ns <- nrow(dat)/2
    iX <- 1:ns
    iY <- (ns+1):(2*ns)
    X <- dat[iX,1,drop=FALSE]
    Y <- dat[iY,1,drop=FALSE]
    rY <- Y - a - b * X
    AA <- omega[iX,iX] + 2*b*omega[iX,iY] + omega[iY,iY]*b^2
    BB <- (omega[iX,iY]+b*omega[iY,iY])%*%rY
    rx <- - solve(AA,BB)
    x <- X - rx
    ry <- dat[iY,1] - a - b*x
    v <- matrix(c(rx,ry),nrow=1,ncol=2*ns)
    (v %*% omega %*% t(v))
}
